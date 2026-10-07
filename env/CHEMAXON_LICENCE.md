# ChemAxon (Certara) academic licence — evaluation

Part of #20. Research only — **no
account was created, no application was submitted, no credentials or personal data were entered.**
This is a report of publicly readable terms plus the exact human action that would be needed next.

Context: #28 fixed this effort as
**internal development / educational use, not commercial distribution.**

## Vendor status

Certara completed its acquisition of Chemaxon on **2 Oct 2024**. All `chemaxon.com` licensing URLs
301-redirect to `certara.com/cxn/...` pages, confirmed live as of this research (14 Aug 2026):

- `chemaxon.com/academic-program` → `https://www.certara.com/cxn/academic-program`
- `chemaxon.com/research-license` → `https://www.certara.com/cxn/research-license`
- `chemaxon.com/academic-offer-form` → `https://www.certara.com/cxn/academic-offer-form`
- `chemaxon.com/research-license-contribution` → `https://www.certara.com/cxn/research-license-contribution`

## What the programme covers

The current Certara pages state the academic programme is for **"academic institutions and the next
generation of scientists"** doing **"non-commercial chemical research"**, offered as:

1. **Marvin Cloud** individual plans (web-based structure editor), with a 14-day free trial and
   "affordable plans for academic users" — pricing not disclosed on-page.
2. An **institutional academic discount** on "a wide range of Chemaxon software," gated on
   "eligible institutions," verified by "talk[ing] to one of our team members."

Historical (pre-Certara-redesign) terms, found via cached/indexed copies of the same URLs and a
ChemAxon forum thread, describe:

- A **free 1-year Individual Research License** for one named individual (thesis/research project),
  non-transferable, contingent on "contribution" (citation in a publication, or acceptable
  alternatives such as a social-media post or short video about the research) tracked over an
  18-month window; failure to contribute blocks renewal.
- Alternatively, a **free-of-charge 2-year academic package**, conditional on crediting ChemAxon in
  any resulting presentation/publication, extendable if a qualifying publication is demonstrated.
- A separate, paid **Commercial Academic License** (~$500/user, historical figure) for
  research groups with a commercial affiliation.

**These historical terms could not be confirmed as still current** — the live Certara pages do not
restate them, and Certara's rewrite may have superseded them. Treat the 1-year / 2-year / $500
figures as **unverified legacy data points**, not current commitments.

## The decisive open question: is the pKa/microspeciation module even in the academic bundle?

ChemAxon's Calculator Plugins are licensed **per module**, not as a single bundle. The current
docs (`docs.chemaxon.com/latest/calculators_licensing.html`) list 7 licensable calculation
groups; the one this effort needs — **pKa, Major Microspecies, and Isoelectric Point** — falls
under the **"Protonation"** license specifically. Free-of-charge calculations are limited to
elemental analysis and 2D polar surface area only. Since Marvin 21.1, in-app (interactive)
property calculations in MarvinSketch no longer require a separate Calculator Plugins licence
beyond the base MarvinBeans licence — **but this exemption explicitly excludes batch/CLI use
(`cxcalc`)**, which is the mode #7 identified as necessary (`cxcalc msdistr` is the only way to
reach fractional microspecies populations; the Python API only exposes `major_microspecies`).

**Neither the current `certara.com/cxn/academic-program` page nor the academic offer form page
names the Protonation module (or the pKa plugin) as included in the free academic bundle.** This
was #7's flagged unknown and it remains unresolved by this research — it is not stated either way
on the pages that are publicly readable without an account.

## Eligibility for an unaffiliated individual

Not established. The gating language is "eligible institutions" / institutional affiliation is
implied throughout; no page states whether a genuinely unaffiliated individual researcher (no
university or company e-mail) qualifies for either the Individual Research License or the
institutional discount. #7 flagged this as open; it is still open after this pass.

## Lead time

**Not published anywhere in the current Certara materials.** Every path — Marvin Cloud individual
plans, the institutional discount, and the academic offer form — routes through a **contact-us /
sales-inquiry form** ("Leave your contact information — Summarize what you are looking for — Our
colleagues will get back to you soon"), with no committed SLA or turnaround time stated. This
confirms #7's finding; no new lead-time data was found.

## What a human would need to do (exact steps, not performed)

1. Go to **https://www.certara.com/cxn/academic-offer-form** (redirect target of
   `chemaxon.com/academic-offer-form`).
2. Fill in institutional/contact details and a description of the research use.
3. Before submitting, get written confirmation (email, not just the form) that the **Protonation /
   pKa Calculator Plugin** specifically — not just base Marvin/MarvinSketch — is included, since
   that module is billed separately in the commercial price list and is not named in the public
   academic-bundle description.
4. Separately ask whether an unaffiliated individual researcher qualifies, since eligibility
   language throughout is institution-centric.
5. Note that #28 already establishes this effort as internal/educational, not commercial
   distribution — that framing should be stated explicitly in the inquiry, since the form's copy
   distinguishes "non-commercial research" from paid commercial/organisational licensing.
6. Expect a sales-team reply with no published SLA; budget unknown lead time (this research could
   not establish a number).

**Recommendation carried forward, matching #7:** do not block on this. The open path (Uni-pKa,
already installed and smoke-tested against 4 compounds in this same ticket — see
`env/SMOKE_TEST_RESULTS.md`) stands on its own and is what's currently wired up. If/when a human
completes steps 1–4 above and confirms the Protonation module is included and lead time, that
becomes the commercial/high-accuracy backend #7 designed the interface to be pluggable with.
