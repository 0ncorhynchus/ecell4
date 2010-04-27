#!/usr/env python

import weakref

import math

import numpy

from utils import *
from surface import *

from gfrdbase import *
import _gfrd

import logging
import myrandom

import itertools

log = logging.getLogger('ecell')

DEFAULT_DT_FACTOR = 1e-5

def calculate_bd_dt(species_list):
    D_list = []
    radius_list = []
    D_max = 0.
    radius_min = numpy.inf
    for species in species_list:
        if D_max < species.D:
            D_max = species.D
        if radius_min > species.radius:
            radius_min = species.radius
    return (radius_min * 2) ** 2 / (D_max * 2)

class BDSimulatorCoreBase(object):
    '''
    BDSimulatorCore borrows the following from the main simulator:
    - species_list
    - reaction_types list (both 1 and 2)
    
    '''
    def __init__(self, main):
        self.main = weakref.proxy(main)

        self.t = 0.0
        self.dt = 0.0

        self.dt_factor = DEFAULT_DT_FACTOR

        self.step_counter = 0
        self.reaction_events = 0

    def initialize(self):
        self.determine_dt()

    def get_next_time(self):
        return self.t + self.dt

    def stop(self, t):
        # dummy
        self.t = t

    def determine_dt(self):
        self.dt = self.dt_factor * \
               calculate_bd_dt(self.main.world.species)
        if __debug__:
            log.debug('bd dt = %g' % self.dt)

    def step(self):
        self.step_counter += 1

        tx = self.main.world.create_transaction()
        ppg = BDPropagator(tx, self.main.network_rules,
                     myrandom.rng, self.dt, self.main.dissociation_retry_moves,
                     [pid for pid, _ in self.main.world])
        while ppg():
            pass

        self.reaction_events += len(ppg.reactions)
        self.t += self.dt

    def check(self):
        for pp in self.tx:
            assert not self.tx.check_overlap(pp)

class BDSimulatorCore(BDSimulatorCoreBase):
    def __init__(self, main):
        BDSimulatorCoreBase.__init__(self, main)

    def initialize(self):
        BDSimulatorCoreBase.initialize(self)

class BDSimulator(ParticleSimulatorBase):
    def __init__(self, world):
        ParticleSimulatorBase.__init__(self, world)
        self.core = BDSimulatorCore(self)
        self.is_dirty = True

    def t(self):
        return self.core.t

    def sett(self, t):
        self.core.t = t

    t = property(t, sett)

    def get_dt(self):
        return self.core.dt

    def get_step_counter(self):
        return self.core.step_counter

    dt = property(get_dt)
    step_counter = property(get_step_counter)


    def initialize(self):
        self.core.initialize()
        self.is_dirty = False

    def get_next_time(self):
        return self.core.t + self.core.dt

    def reset(self):
        # DUMMY
        self.core.t=0

    def stop(self, t):
        # dummy
        self.core.stop(t)

    def step(self):
        self.reaction_type = None

        if self.is_dirty:
            self.initialize()

        self.core.step()

        if __debug__:
            log.info('%d: t=%g dt=%g, reactions=%d, rejected_moves=%d' %
                 (self.step_counter, self.t, self.dt, self.reaction_events,
                  self.rejected_moves))

    def check(self):
        pass
