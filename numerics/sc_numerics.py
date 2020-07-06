#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@project: Stochastic colonization of microbe-free hosts (source code - numerics)
@author: Román Zapién-Campos
(MPI for Evolutionary Biology - zapien@evolbio.mpg.de)
"""


# Import the relevant packages

import os

import numpy as np

import matplotlib.pyplot as mp

from scipy.sparse import lil_matrix

from scipy.sparse.linalg import eigs

from datetime import datetime

from multiprocessing import Pool

from par_numerics import *


# Define working parameters

# Number of cores for parallelization

n_cores = 24

# Numeric precision to achieve

num_prec = 1E-200 #14

# Number of points within a magnitude order

interval_m = 5

interval_t = 5

interval_p = 10

# Ranges of parameters (for migration, host death and frequency in the pool of colonizers)

M, T, P = [], [], []

for i in np.arange(M_min, M_max + 1): M.extend(np.arange(1, 10, 10 / interval_m) * np.power(10., i))

for i in np.arange(T_min, T_max + 1): T.extend(np.arange(1, 10, 10 / interval_t) * np.power(10., i))

for i in np.arange(P_min, P_max + 1): P.extend(np.arange(1, 10, 10 / interval_p) * np.power(10., i))

# Range of states

n = np.arange(0, int(N + 1))


# Build transition probabilities matrix

def Tp_build(N, p, q, m, t, a):

    # Create an empty (sparse) matrix
    
    Tp = lil_matrix((int(N + 1), int(N + 1)))
    
    # Fill the tridiagonal elements
    
    tp_down = (1. - t) * (n / N) * (m * (1. - p) + (1. - m) * ((N - n) / ((1. + a) * (n - 1.) + (N - n))))
    
    tp_up = (1. - t) * ((N - n) / N) * (m * p + (1. - m) * ((1. + a) * n / ((1. + a) * n + (N - n - 1.))))

    tp_remain = 1. - tp_down - tp_up - t
    
    Tp[n[1:], n[1:] - 1] = tp_down[1:]
    
    Tp[n[:-1], n[:-1] + 1] = tp_up[:-1]
    
    Tp[n, n] = tp_remain
    
    # Fill the non tridiagonal elements
        
    # Resetting to state n = 0
    
    if q == 0.:

        Tp[0, 0] = Tp[0, 0] + t

        Tp[1, 0] = Tp[1, 0] + t

        Tp[2:, 0] = t
    
    # Resetting to state n = N
    
    if q == 1.:

        Tp[:-2, -1] = t

        Tp[-2, -1] = Tp[-2, -1] + t 

        Tp[-1, -1] = Tp[-1, -1] + t
    
    # Convert to matrix type CSR and transpose
    
    Tp = Tp.tocsr().T
        
    return Tp


# Build transition rates matrix

def Tr_build(N, p, q, m, t, a):
    
    # Create an empty (sparse) matrix
    
    Tr = lil_matrix((int(N + 1), int(N + 1)))

    # Fill the tridiagonel elements    
    
    tr_down = (1. - t) * (n / N) * (m * (1. - p) + (1. - m) * ((N - n) / ((1. + a) * (n - 1.) + (N - n))))
    
    tr_up = (1. - t) * ((N - n) / N) * (m * p + (1. - m) * ((1. + a) * n / ((1. + a) * n + (N - n - 1.))))
    
    tr_out = - tr_down - tr_up - t
    
    Tr[n[1:], n[1:] - 1] = tr_down[1:]
    
    Tr[n[:-1], n[:-1] + 1] = tr_up[:-1]

    Tr[n, n] = tr_out

    # Fill the non tridiagonal elements
        
    # Resetting to state n = 0
    
    if q == 0.:
        
        Tr[0, 0] = Tr[0, 0] + t
        
        Tr[1, 0] = Tr[1, 0] + t
        
        Tr[2:, 0] = t
    
    # Resetting to state n = N
    
    if q == 1.:
        
        Tr[-1, -1] = Tr[-1, -1] + t
        
        Tr[-2, -1] = Tr[-2, -1] + t
        
        Tr[:-2, -1] = t

    # Convert to matrix type CSR and transpose

    Tr = Tr.tocsr().T
    
    return Tr


# Find the equilibrium of the transition probabilities matrix (Tp) or transition rates matrix (Tr)
    
def T_solve(Tp, Tr):
    
    if Tp == None:
    
        # Find the eigenvalue == 0 of Tr and its associated eigenvector
        
        try:
            
            eig_vv = eigs(Tr, k = 1, sigma = 0, which = 'LM', maxiter = 1E2, tol = 0) #Using sigma destroys the sparcity!
        
        # Fix for singular matrices, adding a small constant to the main diagonal
        
        except:
            
            eig_vv = eigs(Tr + 1E-10 * np.eye(Tr.shape[0]), k = 1, sigma = 0, which = 'LM', maxiter = 1E2, tol = 0) #Using sigma destroys the sparcity!

        # In practice the eigenvalue will be positive, but small. Others will be negative

    if Tr == None:
        
        # Find the eigenvalue == 1 of Tp and its associated eigenvector

        try:
            
            eig_vv = eigs(Tp, k = 1, sigma = 1, which = 'LR', maxiter = 1E2, tol = 0) #Using sigma destroys the sparcity!

        # Fix for singular matrices, adding a small constant to the main diagonal

        except:
                                                                          
            eig_vv = eigs(Tp + 1E-10 * np.eye(Tp.shape[0]), k = 1, sigma = 1, which = 'LR', maxiter = 1E2, tol = 0) #Using sigma destroys the sparcity!

        # In practice the eigenvalue will be positive, but small. Others will be negative

    if Tp != None and Tr != None:
    
        # Find the eigenvalue == 0 of Tr and its associated eigenvector

        try:
            
            eig_vv = eigs(Tr, k = 1, sigma = 0, which = 'LM', maxiter = 1E2, tol = 0) #Using sigma destroys the sparcity!

        except:

            # Fix for singular matrices, adding a small constant to the main diagonal

            try:
                
                eig_vv = eigs(Tr + 1E-10 * np.eye(Tr.shape[0]), k = 1, sigma = 0, which = 'LM', maxiter = 1E2, tol = 0) #Using sigma destroys the sparcity!

            except:

                # Find the eigenvalue == 1 of Tp and its associated eigenvector

                try:
                    
                     eig_vv = eigs(Tp, k = 1, sigma = 1, which = 'LR', maxiter = 1E2, tol = 0) #Using sigma destroys the sparcity!

                # Fix for singular matrices, adding a small constant to the main diagonal

                except:
                    
                     eig_vv = eigs(Tp + 1E-10 * np.eye(Tp.shape[0]), k = 1, sigma = 1, which = 'LR', maxiter = 1E2, tol = 0) #Using sigma destroys the sparcity!

    # Extract the eigenvalue

    eig_val = eig_vv[0]

    # Extract the eigenvector

    eig_vec = eig_vv[1]

    return abs(eig_vec), eig_val


# Compute the probability density function (Phi)

def compute_Phi(arg):
    
    # Unpack parameters
    
    directory, N, p, q, m, t, a, save_files, dim_1 = arg[0], arg[1], arg[2], arg[3], arg[4], arg[5], arg[6], arg[7], arg[8]

    print(p, m, t)

    # Build matrices Tp and Tr assuming infinite living hosts

    if t == None:
        
        Tp = Tp_build(N, p, q, m, 0E0, a)
        
        Tr = Tr_build(N, p, q, m, 0E0, a)

    # Build matrices Tp and Tr assuming finite living hosts
        
    if t != None:
        
        Tp = Tp_build(N, p, q, m, t, a)
        
        Tr = Tr_build(N, p, q, m, t, a)

    # Create array to include initial guesses of the pdf

    Phi = np.zeros((int(N+1), 7))
    
    # Single mode at n = 0
    
    Phi[0, 0] = 1.
    
    # Single mode at n = N
    
    Phi[-1, 1] = 1.
    
    # Single mode at n = N/2
    
    Phi[int(N / 2), 2] = 1.
    
    # Two modes at n = 0 and n = N
    
    Phi[0, 3], Phi[-1, 3] = 0.5, 0.5
    
    # Two modes at n = 0 and n = N/2
    
    Phi[0, 4], Phi[int(N / 2), 4] = 0.5, 0.5
    
    # Two modes at n = N/2 and n = N
    
    Phi[int(N / 2), 5], Phi[-1, 5] = 0.5, 0.5
    
    # Solution given by Tr * Phi = 0 or Tp * Phi = Phi
    
    Phi[:, 6] = T_solve(Tp, Tr)[0].T
    
    Phi[:, 6] = Phi[:, 6] / Phi[:, 6].sum()
    
    # Loop to refine the initial guesses of the equilibrium pdf
    
    i, error = 0, np.ones(7)
    
    while sum(error < num_prec) == 0:
        
        Phi_n = Tp.dot(Phi)
        
        Phi_n = Phi_n/Phi_n.sum(0)
        
        error = abs(Phi_n - Phi).sum(0)
        
        Phi = Phi_n
                
        if i == 1E6: break
    
        i += 1
        
    # Extract the solution out of the initial guesses
    
    ind_min_error = np.argmin(error)
    
    # Renormalize (the numerical error could had been amplified)
    
    Phi = Phi[:, ind_min_error] / sum(Phi[:, ind_min_error])
    
    if save_files == True:
        
        if dim_1 == 'm':
            
            if t != None: 
                
                name_phi = '%sPhi/%s_%s.txt'%(directory, np.where(M == m)[0][0], np.where(T == t)[0][0])
                
                name_error = '%serror/%s_%s.txt'%(directory, np.where(M == m)[0][0], np.where(T == t)[0][0])
            
            else:
                
                name_phi = '%sPhi/%s_None.txt'%(directory, np.where(M == m)[0][0])
            
                name_error = '%serror/%s_None.txt'%(directory, np.where(M == m)[0][0])
            
        elif dim_1 == 'p':
            
            name_phi = '%sPhi/%s_%s.txt'%(directory, np.where(P == p)[0][0], np.where(T == t)[0][0])
            
            name_error = '%serror/%s_%s.txt'%(directory, np.where(P == p)[0][0], np.where(T == t)[0][0])
            
        np.savetxt(name_phi, Phi, delimiter = ',', fmt = '%1.6e')
        
        np.savetxt(name_error, np.array(error[ind_min_error]).reshape(1,), fmt = '%1.6e')
    
    elif save_files == False: return Phi, error[ind_min_error]


# Compute m vs t range

def compute_m_t_range(N, p, q, a, M_min, M_max, T_min, T_max):
    
    # Call the desired number of processors for parallelization
    
    pool = Pool(processes = n_cores)
    
    # Get a name for the directory based on the date and time
    
    c_dir = datetime.now().strftime('%d%m%y%H%M%S')

    # Define main directory

    directory = 'output/%s/'%c_dir

    # Make subdirectories to save Phi, error, and figures

    os.makedirs('output/%s/Phi'%c_dir)
    
    os.makedirs('output/%s/fig'%c_dir)
    
    os.makedirs('output/%s/error'%c_dir)
    
    # Save parameters in text files
    
    np.savetxt('%sN.txt'%directory, np.array(N)[None], fmt = '%1.0e')
    
    np.savetxt('%sp.txt'%directory, np.array(p)[None], fmt = '%1.0e')
    
    np.savetxt('%sq.txt'%directory, np.array(q)[None], fmt = '%1.0e')
    
    np.savetxt('%sa.txt'%directory, np.array(a)[None], fmt = '%1.0e')
    
    np.savetxt('%sMmin.txt'%directory, np.array(M_min)[None], fmt = '%i')
    
    np.savetxt('%sMmax.txt'%directory, np.array(M_max)[None], fmt = '%i')
    
    np.savetxt('%sTmin.txt'%directory, np.array(T_min)[None], fmt = '%i')
    
    np.savetxt('%sTmax.txt'%directory, np.array(T_max)[None], fmt = '%i')
  
    # Compute the pdfs assuming infinite living hosts by parallelizing them
    
    pool.map(compute_Phi, [(directory, N, p, q, m, None, a, True, 'm') for m in M], 1)
    
    # Compute the pdfs by parallelizing them
    
    pool.map(compute_Phi, [(directory, N, p, q, m, t, a, True, 'm') for m in M for t in T], 1)
    
    # Create arrays to store the values to be computed
    
    p0 = np.zeros((len(T), len(M)))
    
    p1N = np.zeros((len(T), len(M)))
    
    pN = np.zeros((len(T), len(M)))
    
    pxi = np.zeros((len(T), len(M)))
    
    dp = np.zeros((len(T), len(M)))
    
    modetype = np.zeros((len(T), len(M)))
    
    error = np.zeros((len(T), len(M)))
    
    # Loop over all the stored pdfs to used them

    for m in range(len(M)):
        
        # Load the m pdf (assuming infinite host lifespan)
        
        Phi_t0 = np.loadtxt('%sPhi/%s_None.txt'%(directory, m), delimiter = ',')
        
        for t in range(len(T)):
            
            # Load the m,t pdf
            
            Phi = np.loadtxt('%sPhi/%s_%s.txt'%(directory, m, t), delimiter = ',')

            # Compute P[x_i < 1/N]

            p0[t, m] = Phi[0]
            
            # Compute P[x_i >= 1/N]
            
            p1N[t, m] = Phi[1:].sum()
            
            # Compute P[x_i >= (N-1)/N]
            
            pN[t, m] = Phi[-1]
            
            # Compute the difference between the model with tau = 0 and tau != 0
            
            dp[t, m] = Dp(Phi, Phi_t0)
            
            # Compute the mode type

            modetype[t, m] = Modetype(Phi)
            
            # Error
            
            error[t, m] = np.loadtxt('%serror/%s_%s.txt'%(directory, m, t), delimiter = ',')
            
            # Plot the pdf to confirm absence of numerical errors

            fig, ax = mp.subplots()
            
            ax.semilogy(Phi)
            
            ax.text(0.05, 0.95, r'error = %1.0e'%error[t, m], transform=ax.transAxes, fontsize=14, verticalalignment='top', bbox = dict(boxstyle='round', facecolor='wheat', alpha=0.5))
            
            ax.set_xlabel(r'$n_i$')
            
            ax.set_ylabel(r'$\Phi(n_i)$')
            
            ax.set_title(r'$m$: %1.0e, $\tau$: %1.0e'%(M[m], T[t]))
            
            mp.savefig('%sfig/m_%1.0e_T_%1.0e.png'%(directory, M[m], T[t]), dpi = 200, bbox_inches = 'tight', format = 'png')
            
            mp.clf()
    
    # Save the arrays of values computed
    
    np.savetxt('%sp0.txt'%directory, p0, delimiter = ',', fmt = '%1.5e') 
    
    np.savetxt('%sp1N.txt'%directory, p1N, delimiter = ',', fmt = '%1.5e')
    
    np.savetxt('%spN.txt'%directory, pN, delimiter = ',', fmt = '%1.5e')
    
    np.savetxt('%spxi.txt'%directory, pxi, delimiter = ',', fmt = '%1.5e')
        
    np.savetxt('%sdp.txt'%directory, dp, delimiter = ',', fmt = '%1.5e')
        
    np.savetxt('%smodetype.txt'%directory, modetype, delimiter = ',', fmt = '%i')
    
    np.savetxt('%serror.txt'%directory, error, delimiter = ',', fmt = '%1.5e')
    

# Compute p vs t range

def compute_p_t_range(N, m, q, a, P_min, P_max, T_min, T_max):
    
    # Call the desired number of processors for parallelization
    
    pool = Pool(processes = n_cores)
    
    # Get a name for the directory based on the date and time
    
    c_dir = datetime.now().strftime('%d%m%y%H%M%S')

    # Define main directory
    
    directory = 'output/%s/'%c_dir
    
    # Make subdirectories to save Phi, error, and figures
    
    os.makedirs('output/%s/Phi'%c_dir)
    
    os.makedirs('output/%s/fig'%c_dir)
    
    os.makedirs('output/%s/error'%c_dir)
    
    # Save parameters in text files
    
    np.savetxt('%sN.txt'%directory, np.array(N)[None], fmt = '%1.0e')
    
    np.savetxt('%sm.txt'%directory, np.array(m)[None], fmt = '%1.0e')
    
    np.savetxt('%sq.txt'%directory, np.array(q)[None], fmt = '%1.0e')
    
    np.savetxt('%sa.txt'%directory, np.array(a)[None], fmt = '%1.0e')
    
    np.savetxt('%sPmin.txt'%directory, np.array(P_min)[None], fmt = '%i')
    
    np.savetxt('%sPmax.txt'%directory, np.array(P_max)[None], fmt = '%i')
    
    np.savetxt('%sTmin.txt'%directory, np.array(T_min)[None], fmt = '%i')
    
    np.savetxt('%sTmax.txt'%directory, np.array(T_max)[None], fmt = '%i')
    
    # Compute the pdfs by parallelizing them
    
    pool.map(compute_Phi, [(directory, N, p, q, m, t, a, True, 'p') for p in P for t in T], 1)
    
    # Create an array to store the values of P[x_i  >= 1/N]
    
    p1N = np.zeros((len(T), len(P)))
    
    error = np.zeros((len(T), len(P)))

    # Loop over all the stored pdfs to used them

    for p in range(len(P)):
                
        for t in range(len(T)):
            
            # Load the p,t pdf
            
            Phi = np.loadtxt('%sPhi/%s_%s.txt'%(directory, p, t), delimiter = ',')

            # Compute P[x_i >= 1/N]

            p1N[t, p] = Phi[1:].sum()
            
            # Error
            
            error[t, p] = np.loadtxt('%serror/%s_%s.txt'%(directory, p, t), delimiter = ',')
            
            # Plot the pdf to confirm absence of numerical errors
            
            fig, ax = mp.subplots()
            
            ax.semilogy(Phi)
            
            ax.text(0.05, 0.95, r'error = %1.0e'%error[t, p], transform=ax.transAxes, fontsize=14, verticalalignment='top', bbox = dict(boxstyle='round', facecolor='wheat', alpha=0.5))
            
            ax.set_xlabel(r'$n_i$')
            
            ax.set_ylabel(r'$\Phi(n_i)$')
            
            ax.set_title(r'$p$: %1.0e, $\tau$: %1.0e'%(P[p], T[t]))
            
            mp.savefig('%sfig/p_%1.0e_T_%1.0e.png'%(directory, P[p], T[t]), dpi = 200, bbox_inches = 'tight', format = 'png')
            
            mp.clf()
    
    # Save the array of values of P[x_i  >= 1/N]
    
    np.savetxt('%sp1N.txt'%directory, p1N, delimiter = ',', fmt = '%1.5e')
    
    np.savetxt('%serror.txt'%directory, error, delimiter = ',', fmt = '%1.5e')
    

# Difference between two probability density functions

def Dp(Phi, Phi_t0):

    # Absolute difference between the two pdfs divided by two
    
    dp = sum(abs(Phi - Phi_t0)) / 2.
    
    return dp


# Find the modes and their location of a probability density function

def Modetype(Phi):
    
    # Set any value smaller than a threshold to 0
    
    Phi[Phi < 1E-9] = 0
    
    # Create arrays to store the values to be computed
    
    modes = []
    
    # Find positive and negative slopes of the function
    
    l1 = np.where((Phi[:-1] - Phi[1:]) > 0)[0]
    
    l2 = np.where((Phi[:-1] - Phi[1:]) < 0)[0]
    
    # Identify the boundaries as modes
    
    if 0 in l1: modes.append(0)
    
    if N-1 in l2: modes.append(int(N))
    
    # Identify any intermediate modes
    
    pe = np.intersect1d(l1, l2 + 1)
    
    for k in pe: modes.append(k)
    
    modes = np.unique(np.array(modes))
    
    # Assign an identifier to each possible combination of modes

    # Unimodal cases

    if len(modes) == 1:
        
        if 0 in modes: modetype = 1
        
        elif N in modes: modetype = 2
        
        elif 0 not in modes and N not in modes: modetype = 3
        
        else: modetype = 0

    # Bimodal cases

    elif len(modes) == 2:
        
        if 0 in modes and N in modes: modetype = 4
        
        elif 0 in modes and N not in modes: modetype = 5
        
        elif 0 not in modes and N in modes: modetype = 6
        
        else: modetype = 0

    # Other cases
        
    else: modetype = 0
    
    return modetype