# CLAUDE.md

A pair of minimal FPsim v3.6 demos illustrating what the FPmod-as-`ss.Pregnancy`
port unlocks when composed with other Starsim modules. Each demo is one script,
one figure, one writeup.

**Purpose is illustrative**, not calibrated for policy decisions.

Two analyses share this repo:

- **Analysis 1 — postpartum HIV/FP one-stop shop** (`demo.py`, spec:
  [`docs/postpartum_one_stop_shop_spec.md`](docs/postpartum_one_stop_shop_spec.md),
  writeup: [`docs/analysis_writeup.md`](docs/analysis_writeup.md)). Shows FPmod
  running in the same sim as STIsim HIV; a single 2-mo postnatal contact
  delivers long-acting contraception + long-acting PrEP.
- **Analysis 2 — postpartum birth spacing → fetal health** (`demo_spacing.py`,
  spec: [`docs/postpartum_spacing_fetal_spec.md`](docs/postpartum_spacing_fetal_spec.md),
  writeup: [`docs/analysis_writeup_spacing.md`](docs/analysis_writeup_spacing.md)).
  Same `PostnatalPackage` intervention (FP arm only), composed with starsim's
  `ssl.mnch.FetalHealth` module and a new `SpacingFetalPenalty` connector that
  reads interpregnancy interval at conception and penalises fetal growth /
  timing when the interval is short.

The `PostnatalPackage` in `demo.py` is the shared intervention — it takes
`offer_fp` / `offer_prep` toggles so both analyses use the same class.

## State of play

**Two illustrative analyses shipped.** Analysis 1 (2026-09-28): 5×2 scenario
grid, HIV-calibrated, ~735k short-interval births averted / ~3.1k female HIV
infections averted at Kenya scale. Analysis 2 (2026-09-29): 5×2 scenario grid,
no HIV, ~487k short-interval births averted / ~2.4k LBW / ~5.2k preterm averted
at Kenya scale.

Two upstream one-liners pushed as branches during analysis 2 (both about the
same class of name-lookup bug already documented for stisim/HIV):

- `starsim@fix/fetalhealth-pregnancy-lookup` — FetalHealth used `sim.demographics.pregnancy` name lookup, which fails for FPmod (name `'fp'`).
- `fpsim@fix/fpmod-fire-delivery-callbacks` — FPmod overrode `_post_delivery` without firing `_delivery_callbacks`.

Neither PR has been opened; branches pushed only.

## Intake

**Model.** `fp.Sim(location='kenya', diseases=sti.HIV(...), networks=[sti.StructuredSexual(), ss.MaternalNet()], ...)`. FPsim on `rc3.6-port` (FPmod as `ss.Pregnancy`); STIsim `v1.6.1`.

**Question.** If a postnatal visit offers long-acting contraception and long-acting PrEP together, what do we get vs. offering either alone?

**Data.** No calibration data — HIV `init_prev`/`beta_m2f` hand-tuned to plausible Kenya prevalence and frozen. Postnatal attendance `p_attend` placeholder 0.5 (should be checked against KDHS postnatal care data if time permits).

**Constraints.** 1-day build; SRH read-out is tomorrow morning. Solo.

## Companion project

[`hiv_kenya`](https://github.com/starsimhub/hiv_kenya) — the calibrated
Kenya HIV model. Distinct research track; no shared code. If Kenya HIV
prevalence needs a better anchor for this demo, borrow the numbers from
hiv_kenya's fit rather than pinging their calibration pipeline.

## Environment

- Python at `/home/robyn/miniconda/bin/python`.
- Both fpsim and stisim installed editable from `gf/*` repos on the
  respective branches noted above. Do not bump either without checking
  the demo still runs.
