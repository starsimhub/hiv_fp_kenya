# hiv_fp_kenya

A pair of minimal FPsim v3.6 demos for Kenya, each a single figure showing
what the FPmod-as-`ss.Pregnancy` port unlocks when composed with other
Starsim modules.

- **Postpartum HIV / FP one-stop shop** — FPsim + STIsim HIV; a 2-mo
  postnatal visit delivers long-acting contraception + long-acting PrEP
  together. Script: [`demo.py`](demo.py); spec:
  [`docs/postpartum_one_stop_shop_spec.md`](docs/postpartum_one_stop_shop_spec.md);
  writeup: [`docs/analysis_writeup.md`](docs/analysis_writeup.md).
- **Postpartum birth spacing → fetal health** — FPsim + starsim's
  `ssl.mnch.FetalHealth`; the same postnatal visit averts short-interval
  pregnancies, and short intervals carry LBW / preterm penalties via a
  small `SpacingFetalPenalty` connector. Script:
  [`demo_spacing.py`](demo_spacing.py); spec:
  [`docs/postpartum_spacing_fetal_spec.md`](docs/postpartum_spacing_fetal_spec.md);
  writeup: [`docs/analysis_writeup_spacing.md`](docs/analysis_writeup_spacing.md).

Illustrative, not calibrated for policy decisions.

Sibling repos: [`hiv_kenya`](https://github.com/starsimhub/hiv_kenya) (the
calibrated Kenya HIV model on stisim 1.7.0),
[`hiv_zambia`](https://github.com/starsimhub/hiv_zambia).
