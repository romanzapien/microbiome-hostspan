#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@project: Stochastic colonization of microbe-free hosts (source code - simulation)
@author: Román Zapién-Campos
(MPI for Evolutionary Biology - zapien@evolbio.mpg.de)
"""


# number of processors to use during the parallelization

n_processors = 24


# import the relevant packages

import numpy as np

import random as rd

import os

import shutil as sh

from time import time

from multiprocessing import Pool

from par_simulation import *


# function to initialize a microbe-free host

def initialize_microbiome(n_microbes):

    microbiome = np.zeros(n_microbe_sp+1, dtype = np.int)

    # generates an 'empty' microbiome

    microbiome[0] = n_microbes

    return microbiome


# function to initialize simulation variables

def initialize(job_name):

    # define array of the species (a list 0, ..., n)

    microbial_species = np.arange(0, n_microbe_sp+1)

    # define the pool of colonizers

    colonizers_pool = np.zeros(n_microbe_sp+1, dtype = np.float)

    if p_symmetry == 'random': colonizers_pool[1:] = np.random.geometric(1E-4, n_microbe_sp)

    elif p_symmetry == 'symmetric': colonizers_pool[1:] = 1.

    colonizers_pool = colonizers_pool/colonizers_pool.sum()

    # define the hosts

    hosts = np.array([initialize_microbiome(n_microbes) for i in np.arange(n_hosts)])

    # save directory name

    _dir = job_name

    return colonizers_pool, microbial_species, hosts, _dir


# function to run a single host simulation

def individual_simulation(args):
    
    # unpack the parameters to use

    microbiome, colonizers_pool, microbial_species = args[0], args[1], args[2]

    # create a vector of the non-neutral parameters

    alphas = np.zeros(n_microbe_sp+1, dtype = np.float)

    # empty-space's

    alphas[0] = alpha

    # initialize arrays to store the longitudinal samples

    microbiome_trajectory = np.zeros((int(t_simulation/t_l_sample),n_microbe_sp+1), dtype = np.int)

    microbiome_t = np.copy(microbiome)

    microbiome_trajectory[0,:] = microbiome_t

    # loop to iterate the death-birth-immigration process of microbes

    for t in np.arange(1, t_simulation):

        # microbial death

        p_death =  1.0 * microbiome_t / n_microbes

        microbiome_t[rd.choices(microbial_species, weights = p_death)] -= 1

        # microbial birth

        p_birth_local = (1.0-m) * (1.+alphas) * microbiome_t / ((1.+alphas) * microbiome_t).sum()

        p_birth = m*colonizers_pool + p_birth_local#/sum(p_birth_local)

        microbiome_t[rd.choices(microbial_species, weights = p_birth)] += 1

        # store data if required

        if t % t_l_sample == 0: microbiome_trajectory[int(t/t_l_sample), :] = microbiome_t

    return microbiome_trajectory


# function to run multiple hosts simulations

def population_simulation(colonizers_pool, microbial_species, hosts, _dir):

    # initialize the multiple processors
    
    pool = Pool(processes = n_processors)

    # create new saving directory

    c_dir = _dir + '_m' + '{0:.1e}'.format(m) + '_a' + '{0:.1e}'.format(alpha)

    os.makedirs('output/'+c_dir)

    # save the parameters file

    sh.copy('par_simulation.py', 'output/'+c_dir+'/')
    
    # create a list of parameters sets

    args = [(hosts[i], colonizers_pool, microbial_species) for i in np.arange(n_hosts)]

    # parallelize the simulations

    hosts_trajectory = pool.map(individual_simulation, iterable = args, chunksize = 1)

    # save the result of the simulations

    grouped_trajectories = [np.zeros((int(t_simulation/t_l_sample), n_hosts)) for n in np.arange(n_microbe_sp+1)]
    
    for h in np.arange(n_hosts):

        for n in np.arange(n_microbe_sp+1): grouped_trajectories[n][:, h] = hosts_trajectory[h][:, n]

    np.savez_compressed('output/%s/data.npz' % c_dir, p = colonizers_pool, empty = grouped_trajectories[0], taxa = grouped_trajectories[1:], m = m, alpha = alpha, n_hosts = n_hosts, n_microbes = n_microbes, t_l_sample = t_l_sample, t_simulation = t_simulation, n_microbe_sp = n_microbe_sp)


# function to set and run a simulation of a host population

def execute(job_name):

    t = time()

    # initialize

    colonizers_pool, microbial_species, hosts, _dir = initialize(job_name)

    # start simulations

    population_simulation(colonizers_pool, microbial_species, hosts, _dir)

    print(time()-t)