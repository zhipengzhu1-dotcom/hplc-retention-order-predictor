#!/usr/bin/env python3
"""What does the model's missing scaffold correlation cost?

#37 measured that per-compound descriptor errors are POSITIVELY correlated
within a scaffold family - rho ~ +0.42 (MeOH), +0.54 (ACN) once propagated to
log k - and warned that drawing compounds independently therefore

    "overstate[s] Var(dlog k) for scaffold pairs by up to ~2x (2s^2(1-rho) vs
     2s^2), which is exactly the mechanism #14 flagged: fictitious order
     uncertainty, every P(A before B) dragged toward 0.5."

The thin slice draws compounds independently (its own honesty list says so), so
that overstatement is live in every number the slice produced - including the
variance attribution that the #30 sweep design was derived from.

This measures the cost. It reuses the slice's own sampler wholesale and changes
exactly one thing: the per-compound descriptor draw becomes

    z_i = sqrt(rho) * z_family + sqrt(1 - rho) * z_i_independent

which induces correlation rho between compounds of a family while leaving each
compound's own 4x4 covariance and marginal spread untouched.

Run with env/model-venv/bin/python (needs numpy, scipy, rdkit).
"""
import itertools
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
SLICE = os.path.join(HERE, "..", "..", "prototype", "thin-slice")
sys.path.insert(0, SLICE)
os.chdir(SLICE)          # slice.py resolves data paths relative to itself

import slice as sl  # noqa: E402
from compounds import COMPOUNDS  # noqa: E402

# Scaffold families among the slice's 16. Names follow the chemistry, not the
# indication: what matters is whether a descriptor error on one member predicts
# an error on another.
FAMILIES = {
    "2-arylpropionic acid": ["Ibuprofen", "Naproxen", "Ketoprofen"],
    "benzoic acid":         ["Benzoic acid", "4-Aminobenzoic acid"],
    "aryloxypropanolamine": ["Propranolol", "Atenolol"],
    "arylamine (tert)":     ["Imipramine", "Lidocaine"],
    "simple aromatic":      ["Toluene", "Naphthalene", "Acetophenone",
                             "Benzamide", "Methylparaben"],
}
FAM_OF = {c: f for f, cs in FAMILIES.items() for c in cs}


def run(rho, seed=sl.SEED, s=sl.S):
    """The slice's sampler with a scaffold-correlated descriptor draw."""
    rng = np.random.default_rng(seed)
    grid, betas, sds = sl.load_column()
    R = sl.design_R()
    stats, _, _, lam, _ = sl.load_ensemble_stats()
    names = [c["name"] for c in COMPOUNDS]
    V = {c["name"]: sl.mcgowan_v(c["smiles"]) for c in COMPOUNDS}

    # ---- per-run draws, identical to the slice
    z6 = rng.standard_normal((s, 6))
    dbeta = np.empty((s, len(grid), 6))
    for g in range(len(grid)):
        L = np.linalg.cholesky(np.outer(sds[g], sds[g]) * R + 1e-12 * np.eye(6))
        dbeta[:, g, :] = z6 @ L.T
    acid_shift = rng.normal(+0.30, 0.30, s)
    base_shift = rng.normal(-0.15, 0.20, s)
    D_run = rng.normal(0.0, 0.30, s)

    # ---- THE ONE CHANGE: shared family component in the descriptor draw
    z_fam = {f: rng.standard_normal((s, 4)) for f in FAMILIES}

    desc, pkas, D = {}, {}, {}
    for c in COMPOUNDS:
        nm = c["name"]
        if c["type"] == "refused":
            continue
        mu, C = stats[nm]
        Cinf = np.outer(lam, lam) * C
        L4 = np.linalg.cholesky(Cinf + 1e-10 * np.eye(4))
        z_ind = rng.standard_normal((s, 4))
        fam = FAM_OF.get(nm)
        if fam is not None and rho > 0:
            z = np.sqrt(rho) * z_fam[fam] + np.sqrt(1.0 - rho) * z_ind
        else:
            z = z_ind
        d4 = mu + z @ L4.T
        d4[:, 2] = np.clip(d4[:, 2], 0.0, None)
        desc[nm] = d4

        pk = (np.tile(np.array(c["pka"], float), (s, 1)) if c["pka"]
              else np.zeros((s, 0)))
        if pk.shape[1]:
            pk = pk + rng.normal(0.0, 0.35, pk.shape)
            if c["type"] == "acid":
                pk[:, 0] += acid_shift
            elif c["type"] == "base":
                pk[:, 0] += base_shift
            elif c["type"] == "ampho":
                pk[:, 0] += base_shift
                pk[:, 1] += acid_shift
        pkas[nm] = pk
        D[nm] = np.clip(1.5 + D_run + rng.normal(0.0, 0.40, s), 0.1, None)

    out = {}
    for phi, pH in sl.OPS:
        beta_phi = np.empty((s, 6))
        for j in range(6):
            y = betas[:, j][:, None] + dbeta[:, :, j].T
            beta_phi[:, j] = sl.PchipInterpolator(grid, y, axis=0)(phi)
        logk = np.full((s, len(names)), np.nan)
        for i, c in enumerate(COMPOUNDS):
            nm = c["name"]
            if c["type"] == "refused":
                continue
            d4 = desc[nm]
            Xs = np.column_stack([np.ones(s), d4[:, 0], d4[:, 1], d4[:, 2],
                                  d4[:, 3], np.full(s, V[nm])])
            logk_n = np.einsum("sj,sj->s", Xs, beta_phi)
            f = sl.f_neutral(c["type"], pkas[nm], pH) if c["pka"] else np.ones(s)
            logk[:, i] = np.log10(10.0 ** logk_n * (f + (1 - f) * 10.0 ** (-D[nm])))
        out[(phi, pH)] = logk
    return names, out


def pair_stats(names, logk):
    """Confident-pair fraction and Var(dlog k), split by same/different family."""
    n = len(names)
    live = [i for i in range(n) if not np.isnan(logk[:, i]).all()]
    conf = same = diff = 0
    var_same, var_diff = [], []
    total = 0
    for i, j in itertools.combinations(live, 2):
        d = logk[:, i] - logk[:, j]
        p = float(np.mean(d < 0))
        p = max(p, 1 - p)
        total += 1
        if p > 0.9:
            conf += 1
        fi, fj = FAM_OF.get(names[i]), FAM_OF.get(names[j])
        if fi is not None and fi == fj:
            same += 1
            var_same.append(float(np.var(d)))
        else:
            diff += 1
            var_diff.append(float(np.var(d)))
    return dict(total=total, conf=conf, n_same=same,
                var_same=float(np.mean(var_same)) if var_same else np.nan,
                var_diff=float(np.mean(var_diff)) if var_diff else np.nan)


def main():
    print("Scaffold correlation: what independent draws cost the slice\n")
    print(f"Families: " + ", ".join(f"{f} ({len(c)})" for f, c in FAMILIES.items()))
    n_same = sum(len(c) * (len(c) - 1) // 2 for c in FAMILIES.values())
    print(f"Within-family compound pairs: {n_same} of 105\n")

    rows = {}
    for rho in (0.0, 0.42, 0.54):
        names, out = run(rho)
        rows[rho] = {op: pair_stats(names, lk) for op, lk in out.items()}

    print(f"{'op (phi, pH)':<16}{'rho':>6}{'confident/105':>15}{'%':>7}"
          f"{'Var same-fam':>14}{'Var diff-fam':>14}")
    for op in sl.OPS:
        for rho in (0.0, 0.42, 0.54):
            r = rows[rho][tuple(op)]
            tag = "  (slice)" if rho == 0.0 else ""
            print(f"{str(tuple(op)):<16}{rho:>6.2f}"
                  f"{f'{r[chr(99)+chr(111)+chr(110)+chr(102)]}/{r[chr(116)+chr(111)+chr(116)+chr(97)+chr(108)]}':>15}"
                  f"{100*r['conf']/r['total']:>6.0f}%"
                  f"{r['var_same']:>14.4f}{r['var_diff']:>14.4f}{tag}")
        print()

    print("Reading it:")
    print("  Var same-fam should FALL as rho rises (2s^2(1-rho)); Var diff-fam")
    print("  should be unmoved. If confident-pair count rises, the slice was")
    print("  under-claiming order confidence on exactly the structurally similar")
    print("  pairs a chromatographer most needs separated.")


if __name__ == "__main__":
    main()
