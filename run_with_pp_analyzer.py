"""
One-off rerun that appends a PostpartumHIVIncidence analyzer to each sim
so we can report HIV incidence *among women within 12 months of delivery*
(the endpoint the spec asks for; stisim's own new_infections_postnatal
counts infant infections via BreastfeedingNet, not adult PP incidence).
"""

import os
os.environ.update(
    OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1',
    NUMEXPR_NUM_THREADS='1', MKL_NUM_THREADS='1',
)

import numpy as np
import pandas as pd
import sciris as sc
import starsim as ss

import demo


class PPHIVInc(ss.Analyzer):
    """Adult female HIV incidence within `window_months` of a live birth."""

    def __init__(self, window_months=12, name='pp_hiv_inc', **kwargs):
        super().__init__(name=name, **kwargs)
        self.window_months = window_months

    def init_results(self):
        super().init_results()
        self.define_results(
            ss.Result('n_pp_women', dtype=int),
            ss.Result('new_pp_infections', dtype=int),
        )

    def step(self):
        sim = self.sim
        fp_mod = sim.demographics.fp
        hiv = sim.diseases.hiv
        ti = self.ti
        # Women in the postpartum window: had a live birth in the last window_months.
        earliest = ti - self.window_months
        in_window = (
            (fp_mod.ti_live_birth >= earliest) & (fp_mod.ti_live_birth <= ti)
            & sim.people.female & sim.people.alive
        )
        pp_uids = ss.uids(np.where(in_window)[0])
        self.results.n_pp_women[ti] = len(pp_uids)
        # New HIV infections this timestep among those women.
        newly = pp_uids[hiv.ti_infected[pp_uids] == ti]
        self.results.new_pp_infections[ti] = len(newly)


END_YEAR = 2040
N_SEEDS = 5


def make_sim(scenario, seed):
    return demo.make_sim(scenario, seed, end_year=END_YEAR, extra_analyzers=[PPHIVInc()])


def run():
    sims = []
    for scen in demo.SCENARIOS:
        for seed in range(N_SEEDS):
            sims.append(make_sim(scen, seed))
    return ss.parallel(sims).sims


if __name__ == '__main__':
    T = sc.timer()
    sims = run()
    T.tt('grid with PPHIVInc')

    rows = []
    for sim in sims:
        scen, seed = sim.label.rsplit('_seed', 1)
        yrs = np.asarray([t.year for t in sim.results.timevec])
        m = (yrs >= 2027) & (yrs <= END_YEAR)
        # Full-population aggregate
        hiv = sim.results.hiv
        # PP-window analyzer results
        pp = sim.results.pp_hiv_inc
        rows.append(dict(
            scenario=scen,
            seed=int(seed),
            new_infections_total=int(np.asarray(hiv.new_infections)[m].sum()),
            new_infections_f=int(np.asarray(hiv.new_infections_f)[m].sum()),
            new_pp_infections=int(np.asarray(pp.new_pp_infections)[m].sum()),
            # Woman-months in PP window (denominator for incidence rate)
            pp_woman_months=int(np.asarray(pp.n_pp_women)[m].sum()),
            live_births=int(np.asarray(sim.results.fp.births)[m].sum()),
            short_intervals=int(np.asarray(sim.results.fp.short_intervals)[m].sum()),
        ))
    df = pd.DataFrame(rows)
    df.to_csv('results/outcomes_pp_analyzer.csv', index=False)
    print(df.to_string(index=False))
