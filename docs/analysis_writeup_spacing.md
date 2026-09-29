# Postpartum birth spacing → fetal health, Kenya — short writeup

## Motivation

Kenya's Demographic and Health Survey records that ~26% of live births follow a preceding birth by less than 24 months, and roughly one in eight by less than 18 months. Both categories carry elevated risk to the newborn: pooled meta-analyses (Conde-Agudelo *et al.* 2006, *JAMA*; DaVanzo 2007) put the adjusted odds ratio at ~1.3–1.5 for low birth weight (LBW, <2500 g) and ~1.3–1.6 for preterm birth (<37 weeks) at short (<18 mo) intervals relative to 24–59 mo. WHO recommends an interpregnancy interval of at least 24 months.

The same 2-month postnatal visit that this repo's HIV+FP analysis exploits also delivers Kenya's most-used contraceptive method — the injectable — to a large share of postpartum women (~43% of first-month postpartum initiators choose injectables). If uptake at that visit lengthens the average interpregnancy interval, the causal chain runs on to fetal outcomes: fewer short-interval pregnancies → fewer growth-restricted and preterm babies.

This analysis asks: *what fraction of the LBW and preterm burden that would otherwise be born in short-interval pregnancies could Kenya avert by scaling a postnatal contraception offer?*

## Intervention modelled

The same `PostnatalPackage` as `docs/postpartum_one_stop_shop_spec.md`, restricted to the FP arm:

| Parameter | Value | Source / notes |
|---|---|---|
| Visit timing | 2 months postpartum | Matches Kenya PNC schedules; avoids FPmod's month-1 contraception overwrite |
| Attendance among eligible | 50% | Placeholder |
| Contraceptive uptake among attendees not already on a method | 30% | Placeholder |
| Method offered | Long-acting injectable | Kenya's most-used method (34% of users); modal postpartum choice (43%) |
| PrEP arm | **Disabled** (`offer_prep=False`) | HIV out-of-scope for this analysis |

Intervention starts 2027; sim runs 2000–2035; reporting window 2027–2035; Kenya-scaled 1:1,353 using the sim-alive vs 2020 Kenya population.

## Model

- **FP engine.** FPsim FPmod, Kenya location, DHS-parameterised. All fertility, contraception choice, and postpartum dynamics from the standard Kenya calibration.
- **Fetal outcomes.** Starsim's `ssl.mnch.FetalHealth` module tracks per-pregnancy growth (`weight_percentile`, `growth_restriction`) and timing (`timing_shift`), and classifies each newborn as LBW / VLBW / SGA at delivery. Preterm classification comes from FPmod / `ss.Pregnancy` (`<37w` GA at birth).
- **Spacing → fetal link.** A small custom connector `SpacingFetalPenalty` (`spacing_fetal.py`) registers a conception callback with FetalHealth. At each conception it looks up months since the mother's most recent live birth and applies:

  | Interval | Growth penalty | Timing shift |
  |---|---|---|
  | < 18 mo | 10% weight reduction | ~Lognormal(2 wk, 0.5 wk) earlier delivery |
  | 18–24 mo | 5% weight reduction | ~Lognormal(1 wk, 0.5 wk) earlier delivery |
  | ≥ 24 mo | none | none |

  Effect sizes chosen so that short-interval pregnancies carry roughly the Conde-Agudelo 2006 ORs. Illustrative, not calibrated.

- **Sim scale.** 5 stochastic replicates at 20,000 agents; runtime ~16 s per pair on 10 workers.
- **No HIV, no STIsim** — clean FPsim + fetal-health story.

## Results

### Birth spacing (Kenya scale, 2027–2035, seed-mean)

- **Baseline live births:** ~15.4 M.
- **Baseline short-interval births (<24 mo):** ~4.4 M (28.5% of live births — consistent with Kenya DHS's 26%).
- **Short-interval births averted by PP-FP:** **~487,000** (range 399k–610k across seeds; **~11% reduction**). Tight, clean signal at 5 seeds.

### Fetal outcomes (Kenya scale, 2027–2035, seed-mean)

- **Baseline LBW births attributable to short spacing:** ~30,300 (0.2% of live births).
- **LBW births averted:** **~2,400** (range 1,450–3,300; **~8% reduction**).
- **Baseline preterm births attributable to short spacing:** ~76,000 (~0.5% of live births).
- **Preterm births averted:** **~5,200** (range 3,400–8,300; **~7% reduction**).
- **Baseline SVN (small-vulnerable newborn = preterm | LBW | SGA):** ~338,000. **SVN averted:** ~28,000 (~8%).

## Key takeaway

At placeholder uptake (50% attend × 30% method), scaling a 2-month postpartum contraception offer to Kenya would avert **~487,000 short-interval births** (11%), **~2,400 LBW babies**, and **~5,200 preterm babies** over 2027–2035, or equivalently about **~28,000 small-vulnerable newborns**. The birth-spacing signal is tight; the fetal-health signal is directionally consistent across every seed but wider — the demo's fetal outcomes represent the spacing-attributable subset only, so the absolute counts are floors, not totals.

## Caveats

- **These are spacing-attributable counts, not totals.** Real Kenya's LBW rate is roughly 10% of live births, not the 0.2% shown here. The demo only applies penalties for short spacing; it does not model maternal nutrition, malaria, HIV, adolescent pregnancy, or any other contributor to LBW / preterm. What's shown is the fraction of the LBW and preterm burden that is *causally attributable* to short intervals in this model.
- **Penalty magnitudes are illustrative.** The 10%/5%/0% growth restriction and 2 wk / 1 wk / 0 wk timing shifts are chosen so that short-interval pregnancies carry Conde-Agudelo-2006-consistent ORs. They are not tuned to Kenya-specific LBW or preterm data. A Kenya-specific calibration would refine these but wouldn't change the shape of the story.
- **Effect combines "fewer pregnancies" and "better-spaced pregnancies".** PP-FP produces ~900k fewer live births as well as better spacing among the pregnancies that occur. The "LBW averted" and "preterm averted" numbers include both channels. A rate-based decomposition (LBW-per-birth in each arm) would separate them; not done here.
- **Uptake parameters (50% × 30%) are placeholders.** Real program attendance and method uptake could push effects up or down substantially.
- **FPsim doesn't tie neonatal mortality to birth weight or GA.** Downstream mortality/DALY effects of averted LBW/preterm are not in-scope of the demo — an FPmod extension or downstream analysis, not this one.
- **Timing-shift is small vs the GA distribution.** Mean GA at birth in the sim is ~40 wk; 1–2 wk shifts move the preterm-rate needle modestly. If the intent were a stronger preterm signal, the timing-shift means would need to increase.

## Sources

- `results/outcomes_spacing.csv` — 5-seed × 20k-agent grid.
- `results/sims_spacing.obj` — pickled sim objects for re-analysis.
- `figures/pp_spacing_fetal.png` — two-panel figure.
- Conde-Agudelo A *et al.* 2006. Birth spacing and risk of adverse perinatal outcomes: a meta-analysis. *JAMA* 295(15):1809–1823.
- Kenya DHS 2022 (birth intervals; contraceptive method mix).
- WHO Report of a Technical Consultation on Birth Spacing, Geneva 2005 (≥24 mo recommendation).

## Upstream fixes required

Two one-line bugs were fixed upstream to make FPsim + FetalHealth compose:

- `starsim@fix/fetalhealth-pregnancy-lookup` — `FetalHealth` used to look up the pregnancy module by attribute name, which fails when FPmod (name `'fp'`) is the pregnancy module. Fixed by using `sim.get_module(ss.Pregnancy)`. Same class of bug as the stisim/HIV lookup already documented in `docs/postpartum_one_stop_shop_spec.md` §5.
- `fpsim@fix/fpmod-fire-delivery-callbacks` — `FPmod._post_delivery` overrode the base hook but never fired `_delivery_callbacks`, so `FetalHealth.on_delivery` never ran. Fixed by adding the callback loop.

Both branches pushed; PRs not opened.
