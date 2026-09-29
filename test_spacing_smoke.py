"""
Smoke test for the birth-spacing -> fetal health analysis.

Two arms, 2 seeds, 5k agents, 2015-2032. Assertions:
- FetalHealth callbacks fire (n_births > 0, sensible mean birth weight).
- SpacingFetalPenalty causes >0 LBW / preterm births at baseline.
- PostnatalPackage reduces short-interval births.
- PostnatalPackage reduces LBW births (directional).
"""

import os
os.environ.update(OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1',
                  NUMEXPR_NUM_THREADS='1', MKL_NUM_THREADS='1')

import numpy as np
import starsim as ss

from demo_spacing import make_sim, N_AGENTS_SMOKE, END_YEAR_SMOKE


def _run_pair(seed):
    base = make_sim('baseline', seed, n_agents=N_AGENTS_SMOKE, end_year=END_YEAR_SMOKE)
    pp = make_sim('pp_fp', seed, n_agents=N_AGENTS_SMOKE, end_year=END_YEAR_SMOKE)
    ss.parallel([base, pp]).sims
    return base, pp


def _summarise(sim):
    fh = sim.custom['fetal_health']
    fp = sim.demographics.fp
    return dict(
        n_births=fh.results.n_births.sum(),
        n_lbw=fh.results.n_lbw.sum(),
        mean_bw=float(fh.results.mean_birth_weight[fh.results.mean_birth_weight > 0].mean()),
        n_short=fp.results.short_intervals.sum(),
        mean_ga=float(fh.results.mean_ga_at_birth[fh.results.mean_ga_at_birth > 0].mean()),
    )


def main():
    seeds = [0, 1]
    base_stats = []
    pp_stats = []
    for s in seeds:
        b, p = _run_pair(s)
        base_stats.append(_summarise(b))
        pp_stats.append(_summarise(p))

    b_births = np.mean([s['n_births'] for s in base_stats])
    b_lbw    = np.mean([s['n_lbw']    for s in base_stats])
    b_bw     = np.mean([s['mean_bw']  for s in base_stats])
    b_short  = np.mean([s['n_short']  for s in base_stats])
    p_lbw    = np.mean([s['n_lbw']    for s in pp_stats])
    p_short  = np.mean([s['n_short']  for s in pp_stats])

    print(f'baseline: n_births={b_births:.0f}, n_lbw={b_lbw:.0f}, mean_bw={b_bw:.0f}g, n_short={b_short:.0f}')
    print(f'pp_fp   : n_lbw={p_lbw:.0f}, n_short={p_short:.0f}')

    assert b_births > 100, f'expected >100 births, got {b_births}'
    assert 2800 < b_bw < 3600, f'implausible mean BW {b_bw}g'
    assert b_lbw > 0, 'connector never fired (no LBW at baseline)'
    assert p_short <= b_short, f'PP-FP did not reduce short intervals ({p_short} vs {b_short})'
    assert p_lbw <= b_lbw, f'PP-FP did not reduce LBW ({p_lbw} vs {b_lbw})'
    print('smoke OK')


if __name__ == '__main__':
    main()
