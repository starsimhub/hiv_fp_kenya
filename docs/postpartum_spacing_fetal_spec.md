# Minimal FPsim + FetalHealth demo — postpartum birth spacing and fetal outcomes, Kenya

**Purpose.** A single figure showing that FPsim 3.6 (FPmod as `ss.Pregnancy`) drives fetal health outcomes through birth spacing, and that a postpartum contraception offer averts short-interval births and — through them — low birth weight and preterm babies. Illustrative, not calibrated, not for decisions.

Companion to `docs/postpartum_one_stop_shop_spec.md` (HIV + FP one-stop shop). Same repo, same intervention class, different second module and different clinical story.

**Repos.** `gf/fpsim` (branch `rc3.6-port`), `gf/starsim` (main). No STIsim in this demo.

## 1. The question

If a postnatal visit offers long-acting injectable contraception, does the resulting reduction in short interpregnancy intervals translate into fewer low-birth-weight and preterm babies?

Why this one:

- It exercises FPsim's real strength — birth spacing — and connects it to a health outcome that policy audiences care about.
- The full causal chain is in the model: intervention → contraception uptake → conception timing → interpregnancy interval → fetal growth / gestational age → LBW / preterm at delivery.
- The intervention class (`PostnatalPackage`) already exists in `demo.py`. The only new modeling is a small connector that reads interpregnancy interval at conception and applies penalties to `FetalHealth`.

## 2. Scenarios (2 × 5 seeds)

| # | Name | Postnatal visit offers |
|---|------|------------------------|
| 0 | Baseline | nothing extra |
| 1 | PP-FP | injectable contraception (`PostnatalPackage(offer_fp=True, offer_prep=False)`) |

Intervention starts in 2027. Sim runs 2000–2035; outcomes reported for 2027–2035, Kenya-scaled to the 2020 population.

## 3. Outputs (one slide, two panels)

- **Panel A — birth spacing:** short birth intervals (<24 mo) averted vs baseline.
- **Panel B — fetal health:** LBW births averted and preterm births averted, side by side.

Table with same numbers, plus seed ranges.

Not reported: DALYs, neonatal mortality (FPsim doesn't tie it to birth weight/GA), any HIV or STI outcome.

## 4. Model setup

```python
import starsim as ss
import starsim.library as ssl
import fpsim as fp
from demo import PostnatalPackage           # from the sibling demo
from spacing_fetal import SpacingFetalPenalty

def make_sim(scenario, seed):
    intvs = []
    if scenario == 'pp_fp':
        intvs.append(PostnatalPackage(offer_fp=True, offer_prep=False))
    return fp.Sim(
        pars=dict(location='kenya', n_agents=20_000,
                  start_year=2000, end_year=2035, rand_seed=seed),
        custom=ssl.mnch.FetalHealth(),
        connectors=[SpacingFetalPenalty()],
        interventions=intvs,
        label=f'{scenario}_seed{seed}',
    )
```

### New connector: `SpacingFetalPenalty(ss.Connector)`

About 40 lines, lives in `spacing_fetal.py`. Registers a conception callback with `FetalHealth`. At each conception:

1. Restrict to parous women (`ti_live_birth.notnan[uids]`).
2. Compute months since last live birth (`sim.ti - ti_live_birth[parous]`; FPmod runs at 1-month timesteps).
3. Bin and apply:
   - **<18 mo:** growth penalty 10%, timing shift ~ Lognormal(mean=2 wk, std=0.5 wk)
   - **18–24 mo:** growth penalty 5%, timing shift ~ Lognormal(mean=1 wk, std=0.5 wk)
   - **≥24 mo:** nothing

Effect sizes chosen so that at Kenya's spacing distribution the baseline LBW rate is a few percentage points above the no-penalty floor and short-interval births carry roughly Conde-Agudelo 2006 odds ratios (~1.4× LBW, ~1.4× preterm, <18 mo vs ≥24 mo). Numbers are illustrative and picked by inspection of the smoke run, not calibrated to Kenya-specific LBW/preterm targets.

### Reused from `demo.py`

- `PostnatalPackage` — refactored to accept `offer_fp` / `offer_prep` toggles so the same class serves both analyses. `offer_prep=False` skips all HIV wiring.
- `_kenya_scale` — same 1:N scale factor from sim alive-count vs Kenya 2020 population.
- `run_grid`, `plot`, `extract_outcomes` structure — parallel to `demo.py`.

## 5. Upstream fixes needed

**starsim `FetalHealth`** (`starsim/library/mnch/fetal_health.py`) looks up the pregnancy module by attribute name (`sim.demographics.pregnancy`, `sim.people.pregnancy`), which fails when FPmod is used (name `'fp'`). Fixed by caching `self.preg = sim.get_module(ss.Pregnancy)` in `init_pre` and using `self.preg` everywhere else, and by using `sim.get_module(ss.Pregnancy)` in `treat_pregnant.step` / `fetal_infection.step`.

Branch: `fix/fetalhealth-pregnancy-lookup`. Pushed, not PRed.

**fpsim `FPmod._post_delivery`** (`fpsim/fpsim/fpmod.py`) overrides the base `_post_delivery` hook but never fires `_delivery_callbacks`, so `FetalHealth.on_delivery` never runs when FPmod is the pregnancy module. Fixed by adding the callback loop at the end of `FPmod._post_delivery`.

Branch: `fix/fpmod-fire-delivery-callbacks`. Pushed, not PRed.

Both are one-liners and belong upstream — same class of bug as the stisim/HIV `'pregnancy' in self.sim.demographics` name lookup documented in `docs/postpartum_one_stop_shop_spec.md` §5.

## 6. Known limitations (found while building)

- **Effect sizes are illustrative.** The three-bin penalty structure is a stylised version of the spacing literature. It gives the right *direction* and roughly the right ORs, but the LBW/preterm counts should not be quoted as Kenya estimates.
- **Only the interval axis.** Growth restriction from other causes (maternal HIV, nutrition, malaria) is not in the demo. The `FetalHealth` module is disease-agnostic and would compose cleanly with an infection-based connector — that's a next step, not this analysis.
- **Ceiling on penalty.** `FetalHealth.apply_growth_restriction` uses diminishing-returns compounding, so repeated short-interval penalties don't send birth weight to zero. Fine for this demo; watch the design if you extend to multi-source penalties.
- **FPsim GA is coarse.** Delivery is scheduled ~9 months from conception with the pregnancy module's GA distribution overlaid; mean GA at birth in the demo is ~40 wk (plausible). Timing shifts of 1–2 wk are small relative to the GA distribution and produce modest preterm-rate shifts.
- **`ti_live_birth` is per-mother, latest-only.** The interval used is time since the *most recent* live birth, not since the last pregnancy end (including losses). Consistent with FPmod's own `short_intervals` result.

## 7. Acceptance checks (smoke)

Smoke test: 2 seeds × baseline vs PP-FP, 5k agents, 2015–2032. Assertions:

- `n_births` from FetalHealth matches FPmod's `births` (delivery callback fires for every live birth).
- Mean birth weight is 2800–3600 g at baseline (plausible).
- `n_lbw > 0` at baseline (connector actually applies penalties).
- PP-FP reduces short-interval births (directional).
- PP-FP reduces LBW births (directional).

Smoke result: mean BW 3280 g, 255 LBW baseline, PP-FP reduces short intervals by ~1% and LBW by ~0.4% at smoke scale — direction correct, magnitude within seed noise; the full grid at 20k × 5 seeds × 8 years of intervention should give visible effects.

## 8. Files

New:
- `spacing_fetal.py` — `SpacingFetalPenalty` connector.
- `demo_spacing.py` — sim builder, run_grid, extract_outcomes, plot.
- `test_spacing_smoke.py` — smoke test.
- `docs/postpartum_spacing_fetal_spec.md` — this file.
- `docs/analysis_writeup_spacing.md` — post-run interpretation.
- `figures/pp_spacing_fetal.png` — output.
- `results/outcomes_spacing.csv`, `results/sims_spacing.obj` — run outputs.

Modified:
- `demo.py` — `PostnatalPackage` takes `offer_fp` / `offer_prep` toggles; HIV wiring guarded by `offer_prep`.

Upstream (branched, pushed, not PRed):
- `starsim@fix/fetalhealth-pregnancy-lookup`
- `fpsim@fix/fpmod-fire-delivery-callbacks`
