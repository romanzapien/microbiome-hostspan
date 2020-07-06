#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@project: Stochastic colonization of microbe-free hosts (parameters - simulation)
@author: Román Zapién-Campos
(MPI for Evolutionary Biology - zapien@evolbio.mpg.de)
"""


# number of hosts

n_hosts = int(5E2)

# number of microbes per host

n_microbes = int(1E4)

# number of species (1 is the minimum number)

n_microbe_sp = int(2E0)

# is the relative abundance in the source community 'symmetric' or 'random'?

p_symmetry = 'symmetric'

# probability of colonization by immigration

m = 1E-2

# Is colonization of empty-space favoured (-1 <= alpha < 0) or unfavoured (1 >= alpha > 0)?

alpha = -0.9

# for how many time-steps does the simulation run?

t_simulation = int(1E7)

 # how often are longitudinal microbiome samples taken?

t_l_sample = int(1E1)