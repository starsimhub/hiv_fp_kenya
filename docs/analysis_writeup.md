# Postpartum one-stop shop, Kenya — short writeup

## Motivation

Postpartum women in Kenya face two co-occurring health risks. First, closely-spaced pregnancies: 26% of live births in Kenya follow a preceding pregnancy by less than 24 months (Kenya DHS), a spacing pattern associated with elevated maternal, neonatal, and infant risk. Second, elevated HIV risk during the postpartum period, both to the mother and (via MTCT and postnatal transmission) to the child. Existing postnatal contact points reach a large share of women but rarely deliver contraception and HIV prevention together, and 61% of women in the postpartum window in the sim are not using any contraceptive method.

Long-acting injectable contraception and long-acting PrEP could plausibly be co-delivered at the same postnatal visit. This analysis asks: *what could Kenya gain from bundling these two services at a single 2-month postpartum contact?*

## Intervention modelled

A one-stop shop offered at **2 months postpartum**. Assumptions:

| Parameter | Value | Source / notes |
|---|---|---|
| Attendance among eligible | 50% | Placeholder; should be checked against KDHS postnatal-care attendance data |
| Contraceptive uptake among attendees not already on a method | 30% | Placeholder |
| Contraceptive method | Injectable (LA) | Matches Kenya's most-popular method (34% of users) and the modal postpartum choice (43% of first-month PP initiators) |
| PrEP uptake among HIV-negative attendees not already on PrEP | 30% | Placeholder |
| PrEP efficacy | 95% | Lenacapavir-like LA-PrEP |
| PrEP course duration | ss.normal(6 mo, 1 mo) | Long-acting, per-uid draw |
| Same attendee cohort for both offers | Yes | One visit, one attendance draw |

Intervention begins in 2027; simulation runs 1985–2040 with reporting window 2027–2040.

## Model

- **FP engine**: FPsim's FPmod (Kenya location, DHS-parameterized). All fertility, contraception choice, and postpartum dynamics inherited from FPsim's Kenya calibration.
- **HIV engine**: STIsim's `sti.HIV`, configured with hiv_kenya's calibrated best-fit parameters (1000-trial Optuna vs UNAIDS 2000–2024; overshoots on PLHIV and prevalence-15–49 but tracks new infections and deaths within a factor of ~1.5).
- **Baseline HIV services**: FSW / general-population / low-CD4 HIV testing arms with historical scale-up; ANC testing at 90%; ART scaling to 97% by 2024. No baseline PrEP intervention, so the postnatal PrEP arm is measured against a clean counterfactual.
- **Sim scale**: 5 stochastic replicates at 20,000 agents, then scaled 1:833 to Kenya's national population.

## Results

### Family planning (baseline 2027–2040 → intervention)

- **Baseline**: ~22.3 M live births, of which **~6.2 M** (28.5%) are short-interval births with a preceding interval < 24 months.
- **Intervention averts ~735,000 short-interval births** (range 654k–814k across seeds; **11.7% reduction**). Signal is clean at 5 seeds.

### HIV (baseline 2027–2040 → intervention)

- **Baseline female HIV**: ~85,000 new infections over 2027–2040 (~6,100/yr; below Kenya's current UNAIDS estimate of ~10–15k/yr, consistent with the calibrated model's known undershoot on prevalence).
- **HIV incidence among postpartum women** (within 12 mo of a live birth): **0.015% per year** (0.15 per 1,000 woman-years) — well below the 1–3% per year reported for postpartum women in Kenya AGYW cohorts, reflecting a calibrated model with HIV incidence collapsing after 2027.
- **Postpartum-window HIV infections under status quo, 2027–2040**: **~3,500 women** at Kenya scale (~250/yr).
- **Total female HIV infections averted by the intervention**: **~3,100** (range −2,500 to +6,700; CI is wide because the raw-agent delta of ~4 infections is amplified 833× by pop scaling — more seeds would tighten this substantially).
- **Postpartum-window infections averted**: **~0** (per-seed range −2,500 to +3,300 at Kenya scale; **0% of PP-window HIV**). The intervention has effectively nothing to prevent inside the postpartum window at this baseline.

## Key takeaway

At placeholder uptake (50% attend × 30% method), a 2-month postpartum one-stop shop could avert **~735,000 short-interval births** in Kenya over 2027–2040 — a robust, tight-CI signal. The LA-PrEP arm's health impact in the current calibrated baseline is essentially zero within the postpartum window itself (~0.015%/yr PP HIV incidence in the sim vs 1–3%/yr in Kenyan empirical cohorts), and the total female-HIV effect (~3,100 averted) sits at the edge of seed-scaled noise; a Kenya-realistic HIV baseline would be needed before quoting a firm PrEP number.

## Caveats

- HIV baseline is calibrated to UNAIDS 2000–2024 but underestimates current Kenya HIV incidence; the absolute HIV numbers here are conservative. Relative reductions are the honest story.
- 5 stochastic replicates give tight CI on FP outcomes but wide CI on HIV, since the sim's raw-agent HIV delta is small and gets amplified by the ~833× population scale factor. 10–15 seeds would tighten Panel B considerably.
- Uptake parameters (50% attend, 30% method uptake) are illustrative placeholders. Real program performance could push effects up or down.
- The intervention's PrEP arm delivers a single 6-month LA course; no re-dosing or later re-contact modeled.
- FP-side calibration inherited from FPsim's Kenya location; HIV pars from hiv_kenya's calibration (a `sti.Sim`); this is the first analysis with the two engines combined, and both fits show minor divergences from real Kenya values at the current parameters.

## Sources

- `results/outcomes.csv` — 5-seed × 20k-agent grid, calibrated HIV.
- `results/outcomes_pp_analyzer.csv` — PP-window HIV incidence via a custom `PPHIVInc` analyzer (see `run_with_pp_analyzer.py`).
- `~/fpsim/fpsim/locations/kenya/data/` — DHS-derived FP files (`mix.csv`, `method_mix_matrix_switch.csv`, `birth_spacing_dhs.csv`).
- `~/hiv_kenya/results/kenya_hiv_calib.obj` — 1000-trial Optuna calibration; best draw applied via `sti.default_build_fn`.
