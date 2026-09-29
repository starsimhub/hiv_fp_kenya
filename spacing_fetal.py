"""
Birth-spacing -> fetal health connector.

Reads months since the mother's most recent live birth at conception and
applies a growth-restriction + timing-shift penalty when the interval is
short. Effect sizes chosen to give roughly Conde-Agudelo 2006 odds ratios
for LBW/preterm at <18 mo vs >=24 mo. Illustrative, not calibrated.

Requires FPmod (which exposes `ti_live_birth`) and starsim's
`ssl.mnch.FetalHealth` in the sim's `custom`.
"""

import numpy as np
import starsim as ss
import starsim.library as ssl


class SpacingFetalPenalty(ss.Connector):
    """Penalise fetal health at conception based on interpregnancy interval."""

    def __init__(self, pars=None, name='spacing_fetal', **kwargs):
        super().__init__(name=name)
        self.define_pars(
            short_cutoff_mo=18.0,
            medium_cutoff_mo=24.0,
            growth_penalty_short=0.10,
            growth_penalty_medium=0.05,
            timing_shift_short=ss.lognorm_ex(mean=ss.weeks(2.0), std=ss.weeks(0.5)),
            timing_shift_medium=ss.lognorm_ex(mean=ss.weeks(1.0), std=ss.weeks(0.5)),
        )
        self.update_pars(pars, **kwargs)

    def init_pre(self, sim):
        super().init_pre(sim)
        if 'fetal_health' not in sim.custom:
            raise ValueError('SpacingFetalPenalty requires ssl.mnch.FetalHealth() in custom.')
        self.fp = sim.get_module(ss.Pregnancy)
        if self.fp is None or not hasattr(self.fp, 'ti_live_birth'):
            raise ValueError('SpacingFetalPenalty requires FPmod (needs ti_live_birth).')
        self.fh = sim.custom['fetal_health']
        self.fh.add_conception_callback(self._on_conception)

    def _on_conception(self, uids):
        parous = uids[self.fp.ti_live_birth.notnan[uids]]
        if not len(parous):
            return

        months = self.sim.ti - self.fp.ti_live_birth[parous]
        short = parous[months < self.pars.short_cutoff_mo]
        medium = parous[(months >= self.pars.short_cutoff_mo) & (months < self.pars.medium_cutoff_mo)]

        if len(short):
            self.fh.apply_growth_restriction(short, self.pars.growth_penalty_short)
            self.fh.apply_timing_shift(short, self.pars.timing_shift_short.rvs(short))
        if len(medium):
            self.fh.apply_growth_restriction(medium, self.pars.growth_penalty_medium)
            self.fh.apply_timing_shift(medium, self.pars.timing_shift_medium.rvs(medium))

    def step(self):
        pass
