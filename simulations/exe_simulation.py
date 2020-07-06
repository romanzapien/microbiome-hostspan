#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@project: Stochastic colonization of microbe-free hosts (execute - simulation)
@author: Román Zapién-Campos
(MPI for Evolutionary Biology - zapien@evolbio.mpg.de)
"""


# import the relevant packages

from datetime import datetime

from sc_simulation import execute


# create a job name for the current execution

job_name = datetime.now().strftime('%d%m%y%H%M%S')

# execute the job

execute(job_name)