# CLAUDE.md

Minimal FPsim + STIsim demo for a Kenya postpartum "one-stop shop" — a
single figure showing that FPsim 3.6 (FPmod as `ss.Pregnancy`) runs in
the same sim as STIsim HIV, with one postnatal intervention that
delivers long-acting contraception + long-acting PrEP together.

**Purpose is illustrative**, not calibrated for HIV and not for decisions.

Full spec — question, scenarios, model setup, the `PostnatalPackage`
intervention, known gotchas (including the stisim `hiv.py` name-lookup
blocker) and the day-plan — is in
[`docs/postpartum_one_stop_shop_spec.md`](docs/postpartum_one_stop_shop_spec.md).
Read that before writing code.

## State of play

**Bootstrap.** Repo scaffolded 2026-09-28. No source code yet. Deliverable is
tomorrow morning's SRH read-out: one two-panel figure comparing baseline
vs PP-FP vs PP-PrEP vs PP-both.

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
