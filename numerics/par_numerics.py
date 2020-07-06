#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@project: Stochastic colonization of microbe-free hosts (parameters - numerics)
@author: Román Zapién-Campos
(MPI for Evolutionary Biology - zapien@evolbio.mpg.de)
"""


# Number of microbes within a host

N = 1E4

# Relative abundance in the pool of colonizing microbes

p = 0E0

# Is the equation representing microbes (q = 0) or unoccupied-space (q = 1)?

q = 1E0

# Is colonization of empty-space favoured (-1 <= alpha < 0) or unfavoured (1 >= alpha > 0)?

a = 0E0

# Probability of colonization by migration

m = 1E-2

# Probability of a death-birth event of hosts

t = 1E-6


# Range of the 'probability of colonization by migration' to explore (log10(m))

M_min, M_max = -7, -1

# Range of the 'probability of a death-birth event of hosts' to explore (log10(t))

T_min, T_max = -10, -3

# Range of the 'relative abundance in the pool of colonizing microbes' to explore (log10(p))

P_min, P_max = -5, -1