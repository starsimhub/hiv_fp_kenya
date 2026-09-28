"""
FPsim + STIsim postpartum "one-stop shop" demo, Kenya.

Illustrative — HIV is hand-tuned, not calibrated. See
docs/postpartum_one_stop_shop_spec.md for scope, gotchas, and scenario
definitions.
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
import stisim as sti
import fpsim as fp


SCENARIOS = ['baseline', 'pp_fp', 'pp_prep', 'pp_both']
N_SEEDS = 3
N_AGENTS = 20_000
START_YEAR = 2000
END_YEAR = 2035
INTERVENTION_START = 2027


class PostnatalPackage(ss.Intervention):
    """Offer LA injectable contraception + LA PrEP at a postnatal visit."""

    def __init__(self, pars=None, name='pp_package', **kwargs):
        super().__init__(name=name)
        self.define_pars(
            offer_fp=True,
            offer_prep=True,
            visit_month=2,
            start_year=INTERVENTION_START,
            p_attend=ss.bernoulli(p=0.5),
            p_fp_uptake=ss.bernoulli(p=0.3),
            p_prep_uptake=ss.bernoulli(p=0.3),
            prep_eff=0.95,
            dur_prep=ss.normal(ss.months(6), ss.months(1)),
        )
        self.update_pars(pars, **kwargs)

    def init_pre(self, sim):
        super().init_pre(sim)
        self.fp_mod = sim.demographics.fp
        self.hiv = sim.diseases.hiv
        self.cm = sim.connectors.contraception
        self.inj_idx = self.cm.methods['inj'].idx
        self._prep_source = abs(hash(self.name)) % 10_000_000

    def init_results(self):
        super().init_results()
        self.define_results(
            ss.Result('n_attended', dtype=int),
            ss.Result('n_fp_started', dtype=int),
            ss.Result('n_prep_started', dtype=int),
        )

    def step(self):
        # Gate on calendar year at the sim level (not the intervention's own timeline).
        if self.sim.t.now('year') < self.pars.start_year:
            return

        fp_mod = self.fp_mod
        hiv = self.hiv
        # ti_live_birth lives on the sim's timeline; compare against sim.ti, not self.ti.
        sim_ti = self.sim.ti
        target_ti = sim_ti - self.pars.visit_month

        alive = self.sim.people.alive
        mask = (fp_mod.ti_live_birth == target_ti) & alive & ~fp_mod.pregnant
        eligible = ss.uids(np.where(mask)[0])
        attendees = self.pars.p_attend.filter(eligible)
        n_att = len(attendees)
        self.results.n_attended[self.ti] = n_att

        if n_att == 0:
            return

        if self.pars.offer_fp:
            not_on_method = attendees[~fp_mod.on_contra[attendees]]
            starters = self.pars.p_fp_uptake.filter(not_on_method)
            if len(starters):
                fp_mod.on_contra[starters] = True
                fp_mod.method[starters] = self.inj_idx
                fp_mod.ever_used_contra[starters] = 1
                # ti_contra is a sim-timeline state; use sim.ti as the base.
                fp_mod.ti_contra[starters] = sim_ti + self.cm.set_dur_method(starters)
            self.results.n_fp_started[self.ti] = len(starters)

        if self.pars.offer_prep:
            hiv_neg = attendees[~hiv.infected[attendees] & ~hiv.on_prep[attendees]]
            starters = self.pars.p_prep_uptake.filter(hiv_neg)
            if len(starters):
                # start_prep expects a scalar ss.Dur; one draw per visit day.
                # ss.normal(ss.months(6), ss.months(1)).rvs() returns floats
                # in years, so wrap in ss.years to match the API.
                dur = ss.years(self.pars.dur_prep.rvs(1)[0])
                hiv.start_prep(
                    starters,
                    eff=self.pars.prep_eff,
                    dur=dur,
                    source_id=self._prep_source,
                    adh=1.0,
                )
            self.results.n_prep_started[self.ti] = len(starters)


def make_sim(scenario, seed, n_agents=N_AGENTS, start_year=START_YEAR, end_year=END_YEAR):
    hiv = sti.HIV(
        init_prev=ss.bernoulli(0.04),
        beta_m2f=0.015,
    )
    intvs = [
        sti.HIVTest(test_prob_data=0.12, name='hiv_test'),
        sti.ART(coverage={'year': [2000, 2010, 2020, 2035], 'value': [0, 0.3, 0.7, 0.9]}),
    ]
    if scenario != 'baseline':
        intvs.append(PostnatalPackage(
            offer_fp=scenario in ('pp_fp', 'pp_both'),
            offer_prep=scenario in ('pp_prep', 'pp_both'),
        ))

    pars = dict(
        location='kenya',
        n_agents=n_agents,
        start_year=start_year,
        end_year=end_year,
        rand_seed=seed,
    )
    return fp.Sim(
        pars=pars,
        diseases=hiv,
        networks=[sti.StructuredSexual(), ss.MaternalNet()],
        interventions=intvs,
        label=f'{scenario}_seed{seed}',
    )


def run_grid(scenarios=SCENARIOS, n_seeds=N_SEEDS, **make_sim_kwargs):
    sims = []
    for scen in scenarios:
        for seed in range(n_seeds):
            sims.append(make_sim(scen, seed, **make_sim_kwargs))
    return ss.parallel(sims).sims


def extract_outcomes(sims, report_start=INTERVENTION_START, report_end=END_YEAR):
    """Per-sim scalar outcomes across the reporting window."""
    rows = []
    for sim in sims:
        scen, seed = sim.label.rsplit('_seed', 1)
        yrs = np.asarray([t.year for t in sim.results.timevec])
        m = (yrs >= report_start) & (yrs <= report_end)
        rows.append(dict(
            scenario=scen,
            seed=int(seed),
            live_births=int(np.asarray(sim.results.fp.births)[m].sum()),
            short_intervals=int(np.asarray(sim.results.fp.short_intervals)[m].sum()),
            new_infections_f=int(np.asarray(sim.results.hiv.new_infections_f)[m].sum()),
            new_infections_postnatal=int(np.asarray(sim.results.hiv.new_infections_postnatal)[m].sum()),
        ))
    return pd.DataFrame(rows)


def plot(df, out='figures/pp_one_stop_shop.png'):
    """Two-panel figure: FP outcomes (Panel A) + HIV outcomes (Panel B)."""
    import pylab as pl
    os.makedirs(os.path.dirname(out), exist_ok=True)

    scen_labels = {
        'baseline': 'Baseline',
        'pp_fp':    'PP-FP',
        'pp_prep':  'PP-PrEP',
        'pp_both':  'PP-both',
    }
    order = [s for s in SCENARIOS if s in df.scenario.unique()]
    x = np.arange(len(order))

    base = df[df.scenario == 'baseline'].mean(numeric_only=True)
    agg = df.groupby('scenario').agg(['mean', 'min', 'max'])

    fig, axes = pl.subplots(1, 2, figsize=(12, 5))

    # Panel A: FP
    ax = axes[0]
    means = [agg.loc[s, ('short_intervals', 'mean')] for s in order]
    lo = [agg.loc[s, ('short_intervals', 'min')] for s in order]
    hi = [agg.loc[s, ('short_intervals', 'max')] for s in order]
    ax.bar(x, means, yerr=[np.array(means) - lo, np.array(hi) - means], capsize=4, color='#4C72B0')
    ax.set_xticks(x)
    ax.set_xticklabels([scen_labels[s] for s in order])
    ax.set_ylabel(f'Short birth intervals <24mo\n({INTERVENTION_START}–{END_YEAR}, cumulative)')
    ax.set_title('Panel A — Family Planning')
    ax.set_ylim(bottom=0)

    # Panel B: HIV
    ax = axes[1]
    means = [agg.loc[s, ('new_infections_f', 'mean')] for s in order]
    lo = [agg.loc[s, ('new_infections_f', 'min')] for s in order]
    hi = [agg.loc[s, ('new_infections_f', 'max')] for s in order]
    ax.bar(x, means, yerr=[np.array(means) - lo, np.array(hi) - means], capsize=4, color='#C44E52')
    ax.set_xticks(x)
    ax.set_xticklabels([scen_labels[s] for s in order])
    ax.set_ylabel(f'New female HIV infections\n({INTERVENTION_START}–{END_YEAR}, cumulative)')
    ax.set_title('Panel B — HIV')
    ax.set_ylim(bottom=0)

    fig.suptitle('Postpartum one-stop shop, Kenya — illustrative', y=1.02)
    fig.tight_layout()
    fig.savefig(out, dpi=120, bbox_inches='tight')
    print(f'saved {out}')
    return fig


if __name__ == '__main__':
    T = sc.timer()
    sims = run_grid()
    T.toc('scenario grid')

    os.makedirs('results', exist_ok=True)
    sc.saveobj('results/sims.obj', sims)

    df = extract_outcomes(sims)
    df.to_csv('results/outcomes.csv', index=False)
    print(df.to_string(index=False))
    plot(df)
