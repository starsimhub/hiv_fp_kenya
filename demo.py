"""
FPsim + calibrated STIsim postpartum "one-stop shop" demo, Kenya.

- HIV configuration and calibrated parameters imported from hiv_kenya
  (see data/kenya_hiv_calib.obj; best-mismatch draw, row 0).
- FP engine is fpsim's FPmod via fp.Sim(location='kenya').
- Two arms: `baseline` (calibrated HIV + hiv_kenya's testing/ART/PrEP
  scale-up) vs `pp_shop` (same, plus a one-stop shop at 2 months
  postpartum that offers LA injectable contraception + LA PrEP to the
  same attendees).

The hiv_kenya calibration overshoots on PLHIV and prevalence but tracks
annual new infections and HIV deaths — appropriate for an illustrative
demo, not a national projection.
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

from interventions import make_hiv_intvs


SCENARIOS = ['baseline', 'pp_shop']
N_SEEDS = 5
N_AGENTS = 20_000
START_YEAR = 1985  # matches hiv_kenya calibration; init_prev_hiv.csv is a 1985-appropriate prevalence
END_YEAR = 2040
INTERVENTION_START = 2027
CALIB_DRAW_IDX = 0


class PostnatalPackage(ss.Intervention):
    """Offer LA injectable contraception and/or LA PrEP at a 2-mo postnatal visit."""

    def __init__(self, pars=None, name='pp_shop', offer_fp=True, offer_prep=True, **kwargs):
        super().__init__(name=name)
        self.offer_fp = offer_fp
        self.offer_prep = offer_prep
        self.define_pars(
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
        self.cm = sim.connectors.contraception
        self.inj_idx = self.cm.methods['inj'].idx
        self.hiv = sim.diseases.get('hiv') if self.offer_prep else None
        if self.offer_prep and self.hiv is None:
            raise ValueError('PostnatalPackage(offer_prep=True) needs sti.HIV in the sim.')
        self._prep_source = abs(hash(self.name)) % 10_000_000

    def init_results(self):
        super().init_results()
        self.define_results(
            ss.Result('n_attended', dtype=int),
            ss.Result('n_fp_started', dtype=int),
            ss.Result('n_prep_started', dtype=int),
        )

    def step(self):
        if self.sim.t.now('year') < self.pars.start_year:
            return

        fp_mod = self.fp_mod
        sim_ti = self.sim.ti
        target_ti = sim_ti - self.pars.visit_month

        eligible = ((fp_mod.ti_live_birth == target_ti) & ~fp_mod.pregnant).uids
        attendees = self.pars.p_attend.filter(eligible)
        self.results.n_attended[self.ti] = len(attendees)

        if len(attendees) == 0:
            return

        if self.offer_fp:
            not_on_method = attendees[~fp_mod.on_contra[attendees]]
            fp_starters = self.pars.p_fp_uptake.filter(not_on_method)
            if len(fp_starters):
                fp_mod.on_contra[fp_starters] = True
                fp_mod.method[fp_starters] = self.inj_idx
                fp_mod.ever_used_contra[fp_starters] = 1
                fp_mod.ti_contra[fp_starters] = sim_ti + self.cm.set_dur_method(fp_starters)
            self.results.n_fp_started[self.ti] = len(fp_starters)

        if self.offer_prep:
            hiv = self.hiv
            hiv_neg = attendees[~hiv.infected[attendees] & ~hiv.on_prep[attendees]]
            prep_starters = self.pars.p_prep_uptake.filter(hiv_neg)
            if len(prep_starters):
                # Per-uid duration draw; wrap in ss.years so the array is Dur-typed
                # (raw numpy floats give a `freq` object when divided by dt).
                dur = ss.years(self.pars.dur_prep.rvs(prep_starters))
                hiv.start_prep(prep_starters, eff=self.pars.prep_eff, dur=dur,
                              source_id=self._prep_source, adh=1.0)
            self.results.n_prep_started[self.ti] = len(prep_starters)


def _load_calib_pars(idx=CALIB_DRAW_IDX):
    calib = sc.loadobj('data/kenya_hiv_calib.obj')
    return calib.df.iloc[idx].to_dict()


def make_sim(scenario, seed, n_agents=N_AGENTS, start_year=START_YEAR, end_year=END_YEAR,
             calib_idx=CALIB_DRAW_IDX, extra_analyzers=None):
    hiv = sti.HIV(
        beta_m2f=0.012,          # overwritten by calib
        eff_condom=0.5,           # overwritten by calib
        init_prev_data=pd.read_csv('data/init_prev_hiv.csv'),
        rel_init_prev=0.5,
        age_bins=[0, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 100],
    )
    nw = sti.StructuredSexual(
        prop_f0=0.79,             # overwritten by calib
        prop_m0=0.75,             # overwritten by calib
        f1_conc=0.15,             # overwritten by calib
        m1_conc=0.15,             # overwritten by calib
        p_pair_form=0.5,          # overwritten by calib
        condom_data=pd.read_csv('data/condom_use.csv'),
    )

    intvs = make_hiv_intvs()
    if scenario == 'pp_shop':
        intvs.append(PostnatalPackage())

    pars = dict(
        location='kenya',
        n_agents=n_agents,
        start_year=start_year,
        end_year=end_year,
        rand_seed=seed,
    )
    sim = fp.Sim(
        pars=pars,
        diseases=hiv,
        networks=[nw, ss.MaternalNet()],
        interventions=intvs,
        analyzers=list(extra_analyzers) if extra_analyzers else None,
        label=f'{scenario}_seed{seed}',
    )
    sim = sti.default_build_fn(sim, _load_calib_pars(calib_idx))
    return sim


def run_grid(scenarios=SCENARIOS, n_seeds=N_SEEDS, **make_sim_kwargs):
    sims = [make_sim(s, k, **make_sim_kwargs) for s in scenarios for k in range(n_seeds)]
    return ss.parallel(sims).sims


KENYA_POP_2020 = 55_000_000  # for scaling sim to national counts


def _kenya_scale(sim, ref_year=2020):
    """fp.Sim doesn't set pop_scale, so results are in raw agents. Compute a
    scale factor from the sim's alive population vs Kenya's real population at
    a reference year."""
    tv = np.asarray([t.year for t in sim.results.timevec])
    idx = int(np.where(tv >= ref_year)[0][0])
    n_alive_sim = float(sim.results.n_alive[idx])
    return KENYA_POP_2020 / n_alive_sim


def extract_outcomes(sims, report_start=INTERVENTION_START, report_end=END_YEAR):
    """Per-sim scalar outcomes across the reporting window, scaled to Kenya."""
    rows = []
    for sim in sims:
        scen, seed = sim.label.rsplit('_seed', 1)
        scale = _kenya_scale(sim)
        # Annualize each flow so sums over years are unambiguous
        new_inf_f = sim.results.hiv.new_infections_f.annualize()
        new_inf = sim.results.hiv.new_infections.annualize()
        births = sim.results.fp.births.annualize()
        short_int = sim.results.fp.short_intervals.annualize()
        # timevec for annualized results is integer-year-valued
        tv = np.asarray(new_inf_f.timevec)
        m = (tv >= report_start) & (tv <= report_end)
        rows.append(dict(
            scenario=scen,
            seed=int(seed),
            scale=scale,
            # Sim-scale raw values
            live_births_sim=births[m].sum(),
            short_intervals_sim=short_int[m].sum(),
            new_infections_f_sim=new_inf_f[m].sum(),
            new_infections_sim=new_inf[m].sum(),
            # Kenya-scale
            live_births=births[m].sum() * scale,
            short_intervals=short_int[m].sum() * scale,
            new_infections_f=new_inf_f[m].sum() * scale,
            new_infections=new_inf[m].sum() * scale,
        ))
    return pd.DataFrame(rows)


def plot(df, out='figures/pp_one_stop_shop.png'):
    """Two-panel: short intervals + female HIV infections averted vs baseline,
    reported at Kenya scale (sim × Kenya_pop / sim_pop_2020)."""
    import pylab as pl
    os.makedirs(os.path.dirname(out), exist_ok=True)

    base = df[df.scenario == 'baseline'].set_index('seed')
    d = df[df.scenario == 'pp_shop'].set_index('seed')
    averted = pd.DataFrame({
        'short_intervals_averted': base['short_intervals'] - d['short_intervals'],
        'infections_averted': base['new_infections_f'] - d['new_infections_f'],
    })

    def bar(ax, metric, color, ylabel):
        vals = averted[metric]
        m, lo, hi = vals.mean(), vals.min(), vals.max()
        ax.bar([0], [m], yerr=[[m - lo], [hi - m]], capsize=6, color=color, width=0.5)
        ax.axhline(0, color='k', lw=0.5)
        ax.text(0, m, f'  {m:+,.0f}\n  (range {lo:,.0f} to {hi:,.0f})',
                va='bottom' if m >= 0 else 'top', ha='left', fontsize=10)
        ax.set_xticks([0])
        ax.set_xticklabels(['PP one-stop shop'])
        ax.set_ylabel(ylabel)
        ax.set_xlim(-0.5, 1.5)

    fig, axes = pl.subplots(1, 2, figsize=(11, 5))
    bar(axes[0], 'short_intervals_averted', '#4C72B0',
        f'Short birth intervals (<24 mo) averted\n({INTERVENTION_START}–{END_YEAR}, Kenya scale)')
    axes[0].set_title('Panel A — Family Planning')
    bar(axes[1], 'infections_averted', '#C44E52',
        f'New female HIV infections averted\n({INTERVENTION_START}–{END_YEAR}, Kenya scale)')
    axes[1].set_title('Panel B — HIV')

    base_si = int(base['short_intervals'].mean())
    base_inf = int(base['new_infections_f'].mean())
    scale = float(base['scale'].mean())
    fig.suptitle(
        f'Postpartum one-stop shop, Kenya — calibrated HIV, scaled 1:{scale:,.0f} to national pop\n'
        f'Baseline mean {INTERVENTION_START}–{END_YEAR}: {base_si:,} short-interval births, '
        f'{base_inf:,} new female HIV infections ({N_SEEDS} seeds)',
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
    sc.saveobj('results/sims.obj', sims)

    df = extract_outcomes(sims)
    df.to_csv('results/outcomes.csv', index=False)
    print(df.to_string(index=False))
    plot(df)
