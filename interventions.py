"""
HIV interventions imported from hiv_kenya: FSW / general / low-CD4 / ANC
testing, ART with future coverage, and PrEP scale-up. The ANC testing
eligibility uses `sim.get_module(ss.Pregnancy)` so it works whether the
pregnancy module is `ss.Pregnancy` (stisim default, name 'pregnancy') or
FPsim's `FPmod` (name 'fp'). Everything else is a straight copy.
"""

import numpy as np
import pandas as pd
import starsim as ss
import stisim as sti


def get_testing_products():
    scaleup_years = np.arange(1990, 2021)
    years = np.arange(1990, 2051)
    n_scaleup = len(scaleup_years)
    n_future = len(years) - n_scaleup

    fsw_prob = np.concatenate([np.linspace(0, 0.75, n_scaleup), np.linspace(0.75, 0.85, n_future)])
    gp_prob = np.concatenate([np.linspace(0, 0.1, n_scaleup), np.linspace(0.1, 0.1, n_future)])
    low_cd4_prob = np.concatenate([np.linspace(0, 0.85, n_scaleup), np.linspace(0.85, 0.95, n_future)])

    def fsw_eligibility(sim):
        return sim.networks.structuredsexual.fsw & ~sim.diseases.hiv.diagnosed & ~sim.diseases.hiv.on_art

    def other_eligibility(sim):
        return ~sim.networks.structuredsexual.fsw & ~sim.diseases.hiv.diagnosed & ~sim.diseases.hiv.on_art

    def low_cd4_eligibility(sim):
        return (sim.diseases.hiv.cd4 < 200) & ~sim.diseases.hiv.diagnosed

    fsw_testing = sti.HIVTest(years=years, test_prob_data=fsw_prob, name='fsw_testing',
                              eligibility=fsw_eligibility, label='fsw_testing')
    other_testing = sti.HIVTest(years=years, test_prob_data=gp_prob, name='other_testing',
                                eligibility=other_eligibility, label='other_testing')
    low_cd4_testing = sti.HIVTest(years=years, test_prob_data=low_cd4_prob, name='low_cd4_testing',
                                  eligibility=low_cd4_eligibility, label='low_cd4_testing')

    def anc_eligibility(sim):
        preg = sim.get_module(ss.Pregnancy)
        return preg.tri1_uids[~sim.diseases.hiv.diagnosed[preg.tri1_uids]]

    anc_testing = sti.HIVTest(test_prob_data=0.9, dt_scale=False, name='anc_testing',
                              eligibility=anc_eligibility, label='anc_testing')

    return fsw_testing, other_testing, low_cd4_testing, anc_testing


def make_hiv_intvs(p_art_projected=0.97):
    """
    Note: hiv_kenya's baseline sti.Prep (scaling to 80% by 2025) is intentionally
    dropped here so the PostnatalPackage's incremental PrEP effect can be
    measured cleanly. This diverges from the calibration setup but is the
    right baseline for the pp_shop counterfactual.
    """
    fsw_testing, other_testing, low_cd4_testing, anc_testing = get_testing_products()

    art_cov = pd.read_csv('data/n_art.csv').set_index('year')
    art_cov['p_art'] = np.nan
    art_cov.loc[2024:, 'p_art'] = p_art_projected
    art = sti.ART(coverage=art_cov)

    return [fsw_testing, other_testing, low_cd4_testing, anc_testing, art]
