# Minimal FPsim + STIsim demo — a postpartum "one-stop shop" in Kenya

**Purpose.** A single figure for tomorrow's SRH read-out showing that FPsim 3.6 (FPmod as `ss.Pregnancy`) runs in the same sim as STIsim HIV, with one intervention acting on both. Timebox: 1 day. It's illustrative, not calibrated for HIV, and not for decisions.

**Repos.** `gf/fpsim` (branch `rc3.6-port`), `gf/stisim` (v1.6.1).

## 1. The question

If a postnatal visit offers long-acting contraception and long-acting PrEP together, what do we get compared with offering either one alone?

Why this one:

- It covers both of your ideas: LA injectable + LA PrEP, delivered at the postpartum contact.
- "Postpartum" is defined natively by FPmod (`sim.people.fp.postpartum`, `ti_delivery`). That's exactly the new 3.6 capability, so the demo showcases the PR directly.
- It uses existing interventions and states on both sides. Only one small custom intervention is needed.
- Fallback variant, same code: swap eligibility to AGYW 15–24 for a general LA-injectable + LA-PrEP co-delivery scenario.

## 2. Scenarios (4 × 5 seeds)

| # | Name | Postnatal visit offers |
|---|------|------------------------|
| 0 | Baseline | nothing extra (status quo FP and HIV programs) |
| 1 | PP-FP | injectable contraception |
| 2 | PP-PrEP | LA-PrEP |
| 3 | PP-both | both, to the same attendees (one visit, one attendance draw) |

The intervention starts in 2027. The sim runs 2000–2035 and outcomes are reported for 2027–2035.

## 3. Outputs (one slide, two panels)

- **Panel A — FP outcomes (FPsim):** short birth intervals (<24 mo) and births averted vs. baseline.
- **Panel B — HIV outcomes (STIsim):** new HIV infections among women in the 12 months after delivery, and total female infections averted vs. baseline.

Also: a one-line table of the same numbers with seed ranges.

Don't report maternal deaths averted. In FPsim, maternal mortality doesn't depend on birth spacing, so the spacing → MMR chain isn't in the model. Short intervals are the honest proxy. Say that out loud; it's a nice "next step".

## 4. Model setup

```python
import numpy as np, sciris as sc, starsim as ss, fpsim as fp, stisim as sti

def make_sim(scenario, seed):
    hiv = sti.HIV(init_prev=..., beta_m2f=...)    # tune so adult female prev ≈ Kenya (check UNAIDS)
    intvs = [sti.HIVTest(...), sti.ART(...)]       # keep HIV roughly stable; defaults OK
    intvs += [sti.ANCTest(visit_prob=0.9)]         # optional; uses get_module(ss.Pregnancy), so works with FPmod
    if scenario != 'baseline':
        intvs += [PostnatalPackage(offer_fp=scenario in ('pp_fp','pp_both'),
                                   offer_prep=scenario in ('pp_prep','pp_both'))]
    return fp.Sim(
        location='kenya', n_agents=20_000, start=2000, stop=2035, rand_seed=seed,
        diseases=hiv,
        networks=[sti.StructuredSexual(), ss.MaternalNet()],   # + ss.BreastfeedingNet() if easy
        interventions=intvs,
    )
```

### Custom intervention: `PostnatalPackage(ss.Intervention)`

About 40 lines. Each step:

1. Find women who delivered `visit_month` timesteps ago (default 2): `fp.ti_delivery == ti - visit_month`, alive, not pregnant.
2. Draw attendance once: `p_attend` (placeholder 0.5, to be checked against KDHS postnatal care data).
3. Among attendees:
   - **FP** (if `offer_fp`): of those not already on a method, `p_fp_uptake` (0.3) start injectables, mirroring the end of `fp.change_initiation.step()` (`fpsim/interventions.py` ~L805–820): `on_contra=True`; `method=cm.methods['inj'].idx`; `ever_used_contra=1`; `ti_contra = ti + cm.set_dur_method(uids)`
   - **PrEP** (if `offer_prep`): of HIV-negative women not on PrEP, `p_prep_uptake` (0.3) call `hiv.start_prep(uids, eff=0.95, dur=ss.months(6), source_id=<unique>, adh=1.0)` (lenacapavir-like; illustrative)
   - **Both:** use the same attendee set for both offers. Optional stretch: correlate uptake with a single willingness draw.
4. Results: `n_attended`, `n_fp_started`, `n_prep_started`.

Use `ss.bernoulli` for the draws (not `np.random`) so scenarios stay common-random-number comparable across seeds.

## 5. Known gotchas (found while reading the code)

**Blocker: STIsim looks up pregnancy by name.** `HIV.include_mtct` is `'pregnancy' in self.sim.demographics` (`stisim/diseases/hiv.py` L208, also L496 and L699 `sim.people.pregnancy`). FPmod is named `'fp'`, so MTCT and pregnancy-aware care-seeking are silently switched off. The same pattern appears in syphilis, BV, PregnancyRiskReduction and ANCSyphTest.

**Fix (preferred, ~30 min):** in `hiv.py`, replace the name lookups with `self.sim.get_module(ss.Pregnancy, die=False)`. `ANCTest` already does this. This is also a good talking point: "the remaining interop work is small and mechanical."

**Hack if stuck:** skip MTCT. The demo outcomes (adult female infections) don't need it.

**Contraception overwrite after delivery.** FPmod sets `ti_contra = ti + 1` at delivery, so the contraception module re-chooses at month 1. Visiting at month 2 (and setting `ti_contra` forward) avoids having the choice overwritten. Check the module step order if you move the visit earlier.

**FPsim sexual activity and the STIsim network are independent.** A woman can conceive without a current partner in `StructuredSexual`, and vice versa. That's fine for a demo, and worth one sentence on the slide ("next step: link them").

**HIV isn't calibrated for Kenya.** Hand-tune `init_prev`/`beta_m2f` for plausible prevalence, then freeze the values. Label the figure "illustrative".

**Check the MaternalNet.** Confirm whether `ss.Pregnancy` adds or needs `ss.MaternalNet` itself, and that `fp.Sim` passes `networks=` through without clobbering anything.

## 6. Day plan and cut lines

| Time | Task | Done when |
|------|------|-----------|
| 1 h | Get `fp.Sim(location='kenya')` + `sti.HIV` + networks running; fix gotcha 1 | Runs end-to-end; `hiv.include_mtct` is True |
| 1 h | Tune HIV to plausible prevalence; run baseline with 5 seeds | Stable prevalence; FP outputs match standalone FPsim |
| 2 h | Write `PostnatalPackage` + a smoke test (starts > 0, targets only postpartum women) | Test passes |
| 1 h | Run 4 × 5 sims in parallel (`ss.parallel` / `sc.parallelize`) | Results saved to disk |
| 1 h | Two-panel figure + table | PNG ready for slide |

If time runs short, cut in this order: `BreastfeedingNet`/MTCT → `ANCTest` → correlated uptake → drop to 3 seeds → drop scenario 1 or 2 (keep baseline vs. both).

## 7. Acceptance checks

- Baseline FP indicators (mCPR, ASFR) match a standalone `fp.Sim(location='kenya')` within noise, so adding HIV doesn't break FP.
- Adding FP-only doesn't change HIV much, and PrEP-only doesn't change births much. These are sanity checks on the plumbing.
- PP-both gives both effects from one set of visits. That's the headline.
