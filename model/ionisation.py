"""Ionisation provider wrapper — the real one, not hardcoded pKa.

The thin slice hardcodes literature pKa values and derives `f_neutral` from a
closed-form sigmoid selected by a hand-assigned `type` ("acid" / "base" /
"ampho" / "refused"). That was deliberate for a prototype — it tests the
composition of the stack, not the pKa provider — but it means the model's
dominant variance term for ionisables is a constant somebody typed in.

This wraps the provider #7 chose and #20 stood up (Uni-pKa via the MIT `unipka`
wrapper) and returns what `spec/scenario-ensemble.md` and #7's interface both
require: **per-microspecies populations with their charges**, not a scalar.

Three things fall out that the hardcoded version could not do:

1. **`f_neutral` is read off the microspecies distribution**, so multiprotic and
   zwitterionic compounds are handled by construction rather than by a
   hand-picked formula per `type`.
2. **The refusal test becomes structural and computable.** CONTEXT.md specifies
   refusal as "a structural test answered from the provider's microspecies
   output" — no neutral microspecies at any pH. The slice hand-labels one
   compound `refused`; here it is derived.
3. **Charges travel with the species**, which #7 requires so a silanol
   ion-exchange term can be attached to the charged fraction later without
   reopening this interface.

Run with env/ionisation-venv/bin/python.
"""
from __future__ import annotations

import json
import os
from dataclasses import dataclass

CACHE = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                     ".ionisation-cache.json")
NEUTRAL_POP_FLOOR = 1e-3     # below this a "neutral" species is not populated


@dataclass(frozen=True)
class Microspecies:
    smiles: str
    charge: int
    fraction: float


@dataclass(frozen=True)
class Speciation:
    """What the provider returns at one pH, for one compound."""
    species: tuple[Microspecies, ...]
    ph: float
    status: str = "ok"          # "ok" | "no_ionisable_site" | "declined"

    @property
    def f_neutral(self) -> float:
        return sum(s.fraction for s in self.species if s.charge == 0)

    @property
    def mean_charge(self) -> float:
        return sum(s.charge * s.fraction for s in self.species)


class Provider:
    """Uni-pKa, with an on-disk cache because inference is slow."""

    def __init__(self, use_simple_smarts: bool = True):
        self._u = None
        self._simple = use_simple_smarts
        self._cache = {}
        if os.path.exists(CACHE):
            with open(CACHE) as f:
                self._cache = json.load(f)

    def _engine(self):
        if self._u is None:
            from unipka import UnipKa
            self._u = UnipKa(use_simple_smarts=self._simple)
        return self._u

    def _save(self):
        with open(CACHE, "w") as f:
            json.dump(self._cache, f)

    def speciate(self, smiles: str, ph: float) -> Speciation:
        """Speciate at one pH.

        ⚠ A provider exception is NOT automatically `PROVIDER_DECLINED`.
        Uni-pKa raises `EnumerationError` when it cannot enumerate two charge
        states — which happens for a compound with **no ionisable site in
        range**, e.g. benzamide. Such a compound is perfectly placeable: it is
        neutral everywhere, `f_neutral` = 1. Mapping that throw onto #33's
        `PROVIDER_DECLINED` tier would send ordinary neutral compounds to the
        expected-but-unplaceable pool, which is exactly wrong.

        The two cases are separated here: `no_ionisable_site` (benign, neutral)
        versus `declined` (genuine, the ~1-in-15 rate #7 measured).
        """
        key = f"{smiles}|{ph:.4f}|{int(self._simple)}"
        if key not in self._cache:
            try:
                df = self._engine().get_distribution(smiles, pH=ph)
            except Exception as exc:                       # provider-specific
                # ⚠ "cannot enumerate 2 charge states" is AMBIGUOUS and the
                # ambiguity is dangerous. It fires both for a neutral compound
                # with no ionisable site (benzamide) AND for a PERMANENTLY
                # CHARGED one with no neutral form (benzyltrimethylammonium) -
                # a quaternary ammonium genuinely has only one charge state.
                # Treating both as benign silently converts a permanently
                # charged compound into a neutral one, which is precisely the
                # error #33/#34/#39 built the refuse tier to prevent.
                #
                # Formal charge separates them, and it is a structural property
                # of the molecule - which is what CONTEXT.md requires the
                # refusal test to be.
                if "enumerate microstates across 2 charge states" in str(exc):
                    from rdkit import Chem
                    m = Chem.MolFromSmiles(smiles)
                    q = Chem.GetFormalCharge(m) if m is not None else None
                    status = "no_ionisable_site" if q == 0 else "permanently_charged"
                else:
                    status = "declined"
                benign = status == "no_ionisable_site"
                self._cache[key] = {"status": status}
                self._save()
                return Speciation(
                    species=((Microspecies(smiles, 0, 1.0),) if benign else ()),
                    ph=ph,
                    status=self._cache[key]["status"])
            cols = {c.lower(): c for c in df.columns}
            sm = cols.get("smiles") or list(df.columns)[0]
            ch = cols.get("charge")
            fr = (cols.get("fraction") or cols.get("population")
                  or cols.get("ratio") or list(df.columns)[-1])
            self._cache[key] = [
                [str(r[sm]), int(r[ch]), float(r[fr])]
                for _, r in df.iterrows()
            ]
            self._save()
        entry = self._cache[key]
        if isinstance(entry, dict):
            st = entry["status"]
            return Speciation(
                species=((Microspecies(smiles, 0, 1.0),) if st == "no_ionisable_site"
                         else ()),
                ph=ph, status=st)
        return Speciation(
            species=tuple(Microspecies(s, c, f) for s, c, f in entry), ph=ph)

    def f_neutral(self, smiles: str, ph: float) -> float:
        return self.speciate(smiles, ph).f_neutral

    def refusal(self, smiles: str, ph_scan=None) -> dict:
        """The #33/#34 structural test, computed rather than hand-labelled.

        Returns the evidence, not just a verdict, so a caller can see WHY: a
        compound with no charge-0 microspecies anywhere in the enumeration is
        structurally refused; one whose neutral form exists but is never
        appreciably populated is the 'degraded' tier, not the refuse tier.
        """
        scan = ph_scan if ph_scan is not None else [0.5 * i for i in range(0, 29)]
        best_ph, best_f = None, 0.0
        neutral_ever_enumerated = False
        for ph in scan:
            sp = self.speciate(smiles, ph)
            if any(s.charge == 0 for s in sp.species):
                neutral_ever_enumerated = True
            f = sp.f_neutral
            if f > best_f:
                best_f, best_ph = f, ph
        st0 = self.speciate(smiles, scan[0]).status
        if st0 == "declined":
            tier = "provider_declined"
        elif st0 == "permanently_charged":
            tier = "refuse"          # no neutral microspecies at any pH (#34)
        elif not neutral_ever_enumerated:
            tier = "refuse"
        elif best_f < NEUTRAL_POP_FLOOR:
            tier = "refuse"          # enumerated but never populated
        elif best_f < 0.5:
            tier = "degraded"
        else:
            tier = "in_envelope"
        return dict(tier=tier, max_f_neutral=best_f, at_ph=best_ph,
                    neutral_enumerated=neutral_ever_enumerated)


def crossover_pka(prov: Provider, smiles: str, lo=0.0, hi=14.0, step=0.02):
    """pKa as the pH where two adjacent net-charge populations cross.

    The textbook definition, and the one env/smoke_test.py found more reliable
    than the provider's own getter labels on zwitterions.
    """
    import numpy as np
    phs = np.arange(lo, hi + 1e-9, step)
    charges = []
    for ph in phs:
        sp = prov.speciate(smiles, float(ph))
        agg = {}
        for s in sp.species:
            agg[s.charge] = agg.get(s.charge, 0.0) + s.fraction
        charges.append(agg)
    states = sorted({c for a in charges for c in a})
    out = []
    for a, b in zip(states, states[1:]):
        prev = None
        for ph, agg in zip(phs, charges):
            # d changes sign where the two adjacent charge states cross. Which
            # WAY it crosses depends on which of the pair is favoured at low pH,
            # so test for any sign change - an earlier version only looked for a
            # downward crossing and silently found nothing at all.
            d = agg.get(a, 0.0) - agg.get(b, 0.0)
            if prev is not None and prev != 0 and (prev > 0) != (d > 0):
                if max(agg.get(a, 0.0), agg.get(b, 0.0)) > 0.01:
                    out.append(round(float(ph), 2))
                break
            prev = d
    return sorted(out)
