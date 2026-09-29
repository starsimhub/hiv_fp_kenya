"""
FPsim + FetalHealth postpartum birth-spacing demo, Kenya.

Illustrative: same PostnatalPackage as demo.py (with `offer_prep=False`),
but the outcomes of interest are short-interval births and the fetal
health consequences of those short intervals (low birth weight, preterm
birth). No HIV, no STIsim.

Chain: postnatal injectable -> fewer short interpregnancy intervals ->
fewer LBW and preterm babies (via SpacingFetalPenalty connector +
starsim FetalHealth).
"""

import os
os.environ.update(
    OMP_NUM_THREADS='1',
    OPENBLAS_NUM_THREADS='1',
    NUMEXPR_NUM_THREADS='1',
    MKL_NUM_THREADS='1',
)

import numpy as np
import pandas as pd
import sciris as sc
import starsim as ss
import starsim.library as ssl
import fpsim as fp

from demo import PostnatalPackage, _kenya_scale
from spacing_fetal import SpacingFetalPenalty


SCENARIOS = ['baseline', 'pp_fp']
N_SEEDS = 5
N_AGENTS = 20_000
START_YEAR = 2000
END_YEAR = 2035
INTERVENTION_START = 2027

N_AGENTS_SMOKE = 5_000
END_YEAR_SMOKE = 2032


def make_sim(scenario, seed, n_agents=N_AGENTS, start_year=START_YEAR, end_year=END_YEAR):
    intvs = []
    if scenario == 'pp_fp':
        intvs.append(PostnatalPackage(offer_fp=True, offer_prep=False))

    pars = dict(
        location='kenya',
        n_agents=n_agents,
        start_year=start_year,
        end_year=end_year,
        rand_seed=seed,
    )
    return fp.Sim(
        pars=pars,
        custom=ssl.mnch.FetalHealth(),
        connectors=[SpacingFetalPenalty()],
        interventions=intvs,
        label=f'{scenario}_seed{seed}',
    )


def run_grid(scenarios=SCENARIOS, n_seeds=N_SEEDS, **make_sim_kwargs):
    sims = [make_sim(s, k, **make_sim_kwargs) for s in scenarios for k in range(n_seeds)]
    return ss.parallel(sims).sims


def extract_outcomes(sims, report_start=INTERVENTION_START, report_end=END_YEAR):
    rows = []
    for sim in sims:
        scen, seed = sim.label.rsplit('_seed', 1)
        scale = _kenya_scale(sim)

        fh = sim.custom['fetal_health']
        fp_mod = sim.demographics.fp

        births = fp_mod.results.births.annualize()
        short = fp_mod.results.short_intervals.annualize()
        n_preterm = fp_mod.results.n_preterm.annualize()
        n_lbw = fh.results.n_lbw.annualize()
        n_svn = fh.results.n_svn.annualize()

        tv = np.asarray(births.timevec)
        m = (tv >= report_start) & (tv <= report_end)

        rows.append(dict(
            scenario=scen,
            seed=int(seed),
            scale=scale,
            live_births_sim=births[m].sum(),
            short_intervals_sim=short[m].sum(),
            n_lbw_sim=n_lbw[m].sum(),
            n_preterm_sim=n_preterm[m].sum(),
            n_svn_sim=n_svn[m].sum(),
            live_births=births[m].sum() * scale,
            short_intervals=short[m].sum() * scale,
            n_lbw=n_lbw[m].sum() * scale,
            n_preterm=n_preterm[m].sum() * scale,
            n_svn=n_svn[m].sum() * scale,
        ))
    return pd.DataFrame(rows)


def plot(df, out='figures/pp_spacing_fetal.png'):
    import pylab as pl
    os.makedirs(os.path.dirname(out), exist_ok=True)

    base = df[df.scenario == 'baseline'].set_index('seed')
    d = df[df.scenario == 'pp_fp'].set_index('seed')

    averted = pd.DataFrame({
        'short_intervals_averted': base['short_intervals'] - d['short_intervals'],
        'lbw_averted':             base['n_lbw']           - d['n_lbw'],
        'preterm_averted':         base['n_preterm']       - d['n_preterm'],
    })

    def bar(ax, metric, color, ylabel):
        vals = averted[metric]
        m, lo, hi = vals.mean(), vals.min(), vals.max()
        ax.bar([0], [m], yerr=[[max(0, m - lo)], [max(0, hi - m)]], capsize=6, color=color, width=0.5)
        ax.axhline(0, color='k', lw=0.5)
        ax.text(0, m, f'  {m:+,.0f}\n  (range {lo:,.0f} to {hi:,.0f})',
                va='bottom' if m >= 0 else 'top', ha='left', fontsize=10)
        ax.set_xticks([0])
        ax.set_xticklabels(['PP-FP'])
        ax.set_ylabel(ylabel)
        ax.set_xlim(-0.5, 1.5)

    fig, axes = pl.subplots(1, 2, figsize=(12, 5))
    bar(axes[0], 'short_intervals_averted', '#4C72B0',
        f'Short birth intervals (<24 mo) averted\n({INTERVENTION_START}-{END_YEAR}, Kenya scale)')
    axes[0].set_title('Panel A - Birth spacing')

    ax2 = axes[1]
    lbw_v = averted['lbw_averted']
    pt_v = averted['preterm_averted']
    for i, (name, vals, color) in enumerate([
        ('LBW births\naverted', lbw_v, '#8172B2'),
        ('Preterm births\naverted', pt_v, '#CCB974'),
    ]):
        m, lo, hi = vals.mean(), vals.min(), vals.max()
        ax2.bar([i], [m], yerr=[[max(0, m - lo)], [max(0, hi - m)]], capsize=6, color=color, width=0.5)
        ax2.text(i, m, f'  {m:+,.0f}', va='bottom' if m >= 0 else 'top', ha='center', fontsize=10)
    ax2.axhline(0, color='k', lw=0.5)
    ax2.set_xticks([0, 1])
    ax2.set_xticklabels(['LBW', 'Preterm'])
    ax2.set_ylabel(f'Births averted ({INTERVENTION_START}-{END_YEAR}, Kenya scale)')
    ax2.set_title('Panel B - Fetal health')

    base_si = int(base['short_intervals'].mean())
    base_lbw = int(base['n_lbw'].mean())
    base_pt = int(base['n_preterm'].mean())
    scale = float(base['scale'].mean())
    fig.suptitle(
        f'Postpartum birth spacing -> fetal health, Kenya\n'
        f'Baseline {INTERVENTION_START}-{END_YEAR}: {base_si:,} short-interval births, '
        f'{base_lbw:,} LBW, {base_pt:,} preterm ({N_SEEDS} seeds; scale 1:{scale:,.0f})',
        y=1.02, fontsize=11)
    fig.tight_layout()
    fig.savefig(out, dpi=120, bbox_inches='tight')
    print(f'saved {out}')
    return fig


if __name__ == '__main__':
    T = sc.timer()
    sims = run_grid()
    T.toc('scenario grid')

    os.makedirs('results', exist_ok=True)
    sc.saveobj('results/sims_spacing.obj', sims)

    df = extract_outcomes(sims)
    df.to_csv('results/outcomes_spacing.csv', index=False)
    print(df.to_string(index=False))
    plot(df)
