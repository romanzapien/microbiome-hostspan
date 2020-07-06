#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@project: Stochastic colonization of microbe-free hosts (execute - numerics)
@author: Román Zapién-Campos
(MPI for Evolutionary Biology - zapien@evolbio.mpg.de)
"""


# import the relevant packages

from sc_numerics import *

# compute the numeric solution of the model in a given range

compute_m_t_range(N, p, q, a, M_min, M_max, T_min, T_max)

compute_p_t_range(N, p, q, a, P_min, P_max, T_min, T_max)