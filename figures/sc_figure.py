#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@project: Stochastic colonization of microbe-free hosts (source-code - figures)
@author: Román Zapién-Campos
(MPI for Evolutionary Biology - zapien@evolbio.mpg.de)
"""


# Import the relevant packages

import numpy as np

import mpmath as mpm

from numba import jit

import matplotlib.pyplot as mp

from mpl_toolkits.axes_grid1.inset_locator import inset_axes

import gc


# Sample times from an exponential decaying distribution

@jit(nopython = True)
def time_sample(time_par):

    return np.random.exponential(time_par)


# Figure of the colonization simulations

def fig_simulation(ax, data_dir, title_text, label_text, y_axis, T, T_label):

    # Parameters
    
    n_hosts = 5E2

    n_microbes = 1E4

    t_l_sample = 1E1

    t_simulation = 1E7

    n_microbe_sp = 2E0
    
    start_spared_tsteps = 1E4

    sampled_tsteps = 5E3

    n_sampled_hosts = 5E1

    n_bins = 30

    rep = 1E1
    
    # Load data
    
    data = np.load(data_dir+'data.npz')

    data_files = data['taxa']
    
    # Calculate average trajectory
    
    counts_mean = np.zeros((int(t_simulation/t_l_sample), int(n_microbe_sp)))

    for i in np.arange(n_microbe_sp, dtype = np.int): counts_mean[:, i] = data_files[i].sum(1)/(n_microbes*n_hosts)

    # Set figure axes limits

    ax.set_xlim(-1E-2, 1+1E-2)

    ax.set_ylim(-1E-2, 1+1E-2)
    
    # Sample a subset of hosts randomly
    
    sampled_hosts = np.random.choice(int(n_hosts), int(n_sampled_hosts), replace = False)

    # Plot the colonization trajectory of the hosts sampled

    for j in sampled_hosts:
        
        # Every datapoint before 'start_spared_tsteps' is plotted, but after this only some are plotted
        # this reduces the size of the figure and limits overcrowding
        
        # create random time-points list after 'start_spared_tsteps'
        
        sampled_ts = np.arange(start_spared_tsteps, t_simulation/t_l_sample, dtype = np.int)

        np.random.shuffle(sampled_ts)

        # join to points before 'start_spared_tsteps'

        sampled_ts = np.hstack((np.arange(start_spared_tsteps, dtype = np.int), sampled_ts[:int(sampled_tsteps)]))

        sampled_ts.sort()

        # plot according to the color code

        color = np.arange(t_simulation/t_l_sample)[sampled_ts]

        ax.scatter(data_files[0][sampled_ts, j]/n_microbes, data_files[1][sampled_ts, j]/n_microbes, cmap = mp.cm.gnuplot, vmin = 0, vmax = t_simulation/t_l_sample, s = 0.1, c = color, rasterized = True)

    # Plot the average trajectory

    ax.scatter(counts_mean[:, 0], counts_mean[:, 1], s = 0.1, c = 'green', rasterized = True)

    # Plot the carrying capacity line

    ax.plot(np.linspace(0., 1., 50), np.linspace(1., 0., 50), color = 'gray', linestyle = '--', zorder = 0)

    # Set x label

    ax.set_xlabel(r'microbial taxon 1 ($x_1$)', fontsize = 16)

    # Set y label

    if y_axis == True: ax.set_ylabel(r'microbial taxon 2 ($x_2$)', fontsize = 16)

    if y_axis != True: ax.yaxis.set_ticklabels([])

    # Set title

    ax.set_title(title_text, fontsize = 16, fontweight = 'bold', pad = 10)

    # Set panel label

    ax.text(-0.03, 1.03, label_text, fontsize = 25, fontweight = 'bold', transform = ax.transAxes)
    
    # Create an inset panel
    
    inset = inset_axes(ax, width = '40%', height = '40%', loc = 'upper right', borderpad = 1)
    
    # Sample times of host death according to the host death rate 1/T
    
    host_db_time = []

    for i in np.arange(int(n_hosts*rep)): host_db_time.append(int(time_sample(1./T)/t_l_sample))

    # Create a list of hosts to die according to those times, randomly

    host_id = np.zeros(int(n_hosts*rep), dtype = np.int)

    for r in np.arange(rep, dtype = np.int):

        seq = np.arange(n_hosts)

        np.random.shuffle(seq)

        host_id[int(n_hosts*r):int(n_hosts*(1+r))] = seq
    
    # Subsample the host population as many times as indicated by 'rep'
    
    subsample = np.zeros((int(n_microbe_sp), int(n_hosts*rep)))
    
    for j in np.arange(int(n_microbe_sp), dtype = np.int):

        subsample[j, :] = data_files[j][host_db_time, host_id]
    
    subsample = subsample/(rep*n_microbes)

    # Plot histogram of relative frequency of the focal microbe in the host population

    inset.hist(subsample[0, :]/subsample.sum(0), bins = n_bins, range = (0, 1), color = 'green', edgecolor = 'black', alpha = 1.)

    # Set inset title, axes ticks and labels

    inset.text(0.5, 0.9, r'$\tau = $'+T_label, horizontalalignment = 'center', verticalalignment = 'center', transform = inset.transAxes, fontsize = 12)

    inset.set_xticks([0., 0.5, 1.])

    inset.set_yticklabels(['', '', '', '', '', ''])

    inset.set_ylabel(r'frequency', fontsize = 12)

    inset.set_xlabel(r'$x_1/(x_1 + x_2)$', fontsize = 12)

    # Set inset axes limits

    inset.set_xlim(0, 1)

    inset.set_ylim(0, 0.5*n_hosts*rep)

    return ax


# Figure of the difference between two distributions

def fig_dp(ax, title_text, label_text, y_axis, data, c_levels,c_levels_labels):

    # Load data
    
    Z_dp = np.loadtxt(data+'dp.txt', delimiter = ',')

    M_min, M_max = np.loadtxt(data+'Mmin.txt'), np.loadtxt(data+'Mmax.txt')

    T_min, T_max = np.loadtxt(data+'Tmin.txt'), np.loadtxt(data+'Tmax.txt')

    # Create ranges of migration 'm' and host death 't' rates 

    M, T = [], []

    for i in np.arange(M_min, M_max+1): M.extend(np.arange(1, 10, 2)*np.power(10., i))

    for i in np.arange(T_min, T_max+1): T.extend(np.arange(1, 10, 2)*np.power(10., i))

    # Compute the log10 of the data

    Z_dp = np.log10(np.abs(Z_dp))        

    # Create 2D space

    M_, T_ = np.meshgrid(M, T)       

    # Set log10 scale of axes

    ax.set_xscale('log')

    ax.set_yscale('log')
    
    # Plot the difference btw. distributions

    ax.pcolormesh(M_, T_, Z_dp, cmap = mp.cm.summer, linewidth = 1E-8, zorder = 0)

    # Draw the contour lines and label them

    cax = ax.contour(M_, T_, Z_dp, levels = c_levels, linewidths = 2.0, colors = 'k', origin = 'lower', linestyles = 'solid')

    mp.clabel(cax, inline = 1, fontsize = 10, inline_spacing = 8, use_clabeltext = True, fmt = c_levels_labels, manual = True)

    # Set title, panel label and axes labels

    ax.set_title(title_text, fontsize = 16, fontweight = 'bold', pad = 10)

    ax.text(-0.03, 1.03, label_text, fontsize = 25, fontweight = 'bold', transform = ax.transAxes)

    if y_axis == True: ax.set_ylabel(r'prob. of host death ($\tau$)', fontsize = 16)

    if y_axis != True: ax.yaxis.set_ticklabels([])

    ax.set_xlabel(r'prob. of migration ($m$)', fontsize = 16)
        
    return ax


# Figure of the modality of a distribution

def fig_mode(ax, title_text, label_text, y_axis, data, cmap, norm):

    # Load data
    
    Z_mode = np.loadtxt(data+'modetype.txt', delimiter=',')

    M_min, M_max = np.loadtxt(data+'Mmin.txt'), np.loadtxt(data+'Mmax.txt')

    T_min, T_max = np.loadtxt(data+'Tmin.txt'), np.loadtxt(data+'Tmax.txt')

    # Create ranges of migration 'm' and host death 't' rates 

    M, T = [], []

    for i in np.arange(M_min, M_max+1): M.extend(np.arange(1, 10, 2)*np.power(10., i))

    for i in np.arange(T_min, T_max+1): T.extend(np.arange(1, 10, 2)*np.power(10., i))

    # Plot the modality of distributions

    cax = ax.imshow(Z_mode, origin = 'lower', cmap = cmap, norm = norm, vmin = 1, vmax = 6)

    # Set axes ticks and tick labels using a log10 scale

    ax.set_xticks(np.arange(0, len(M), 5))

    ax.set_yticks(np.arange(0, len(T), 5))

    ax.set_xticklabels([r'$10^{'+str(int(i))+'}$' for i in np.arange(M_min, M_max+1)])

    ax.set_yticklabels([r'$10^{'+str(int(i))+'}$' for i in np.arange(T_min, T_max+1)])

    # Set title, panel label and axes labels

    ax.set_title(title_text, fontsize = 16, fontweight = 'bold', pad = 10)

    ax.text(-0.03, 1.03, label_text, fontsize = 25, fontweight = 'bold', transform = ax.transAxes)

    if y_axis == True: ax.set_ylabel(r'prob. of host death ($\tau$)', fontsize = 16)

    if y_axis != True: ax.yaxis.set_ticklabels([])

    ax.set_xlabel(r'prob. of migration ($m$)', fontsize = 16)

    return ax


# Figure of the probability

def fig_p(ax, title_text, label_text, y_axis, data, c_p, c_levels, c_levels_labels):
    
    # Load data
    
    Z_p = np.loadtxt(data+'%s.txt'%c_p, delimiter = ',')

    M_min, M_max = np.loadtxt(data+'Mmin.txt'), np.loadtxt(data+'Mmax.txt')

    T_min, T_max = np.loadtxt(data+'Tmin.txt'), np.loadtxt(data+'Tmax.txt')

    # Create ranges of migration 'm' and host death 't' rates 

    M, T = [], []

    for i in np.arange(M_min, M_max+1): M.extend(np.arange(1, 10, 2)*np.power(10., i))

    for i in np.arange(T_min, T_max+1): T.extend(np.arange(1, 10, 2)*np.power(10., i))

    # Set the values of probability smaller than a threshold, equal to zero

    Z_p[Z_p < 1E-12] = 1E-323

    # Compute the log10 of the data

    Z_p = np.log10(np.abs(Z_p))

    # Create 2D space

    M_, T_ = np.meshgrid(M, T)        

    # Set log10 scale of axes

    ax.set_xscale('log')

    ax.set_yscale('log')

    # Create shadow to indicate where the probability is greater than 0.5

    Z_p_s = Z_p < mpm.log(mpm.mpf('0.5'), 10.)

    Z_p_s = np.ma.masked_where(Z_p_s == True, Z_p_s)

    # Plot the probability

    ax.pcolormesh(M_, T_, Z_p, cmap = mp.cm.summer, linewidth = 1E-8, zorder = 0, rasterized=True)

    # Shadow the regions where the probability is larger than 0.5

    sax = ax.pcolormesh(M_, T_, Z_p_s, cmap = mp.cm.gray, alpha = 0.12, linewidth = 1E-8, zorder = 1)

    # Draw the contour lines and label them

    cax = ax.contour(M_, T_, Z_p, levels = c_levels, linewidths = 2.0, colors = 'k', origin = 'lower', linestyles = 'solid')

    mp.clabel(cax, inline = 1, fontsize = 10, inline_spacing = 5, use_clabeltext = True, fmt = c_levels_labels, manual = True)

    # Set title, panel label and axes labels

    ax.set_title(title_text, fontsize = 16, fontweight = 'bold', pad = 10)

    ax.text(-0.03, 1.03, label_text, fontsize = 25, fontweight = 'bold', transform = ax.transAxes)

    if y_axis == True: ax.set_ylabel(r'prob. of host death ($\tau$)', fontsize = 16)

    if y_axis != True: ax.yaxis.set_ticklabels([])

    ax.set_xlabel(r'prob. of migration ($m$)', fontsize = 16)

    return ax


# Figure of the simulations vs. numerics

def fig_simvsmodel(ax, title_text, label_text, y_axis, num, sim, c_p, c_t, y_label):

    # Parameters
    
    n_hosts = 5E2

    n_microbes = 1E4

    t_l_sample = 1E1

    t_simulation = 1E7

    n_microbe_sp = 2E0

    M_min, M_max = -7, -1

    T_min, T_max = -10, -3

    interval_m, interval_t = 5
    
    rep = 1E0

    # Create ranges of migration 'm' and host death 't' rates 
    
    M, T = [], []

    for i in np.arange(M_min, M_max + 1): M.extend(np.arange(1, 10, 10 / interval_m) * np.power(10., i))

    for i in np.arange(T_min, T_max + 1): T.extend(np.arange(1, 10, 10 / interval_t) * np.power(10., i))

    M, T = np.array(M), np.array(T)

    # Load simulation data

    sim_m = ['../data/simulations/m1.0e-0%i_%s/'%(i,sim) for i in range(7, 0, -1)]

    # Load numeric-solution data
    
    num_p = np.loadtxt(num+'%s.txt'%c_p, delimiter = ',')

    # Define subspace to plot

    m_ind = np.arange(0, 35, 5, dtype = np.int)

    t_ind = np.arange(20, 36, 1, dtype = np.int)
    
    # Define color of lines
    
    colors = ['brown', 'red', 'orange', 'gold', 'green', 'blue', 'turquoise']
    
    # Differentiate btw. plots of microbial taxa and empty-space
    
    if c_t == 'empty-space': ind = 0

    if c_t == 'taxa': ind = 1
    
    # Span the range of migration rates 'm'
    
    for m in range(len(m_ind)):
        
        # Load specific data
        
        data_dir = sim_m[m]
        
        # Create arrays to store the probabilites of simulation and numeric data
        
        p_sim = np.zeros(len(t_ind))

        p_num = np.zeros(len(t_ind))
        
        # Span the range of host death rates 't'
        
        for t in range(len(t_ind)):
    
            # Draw time for host death
            
            host_db_time = []

            i = 0

            while i < int(n_hosts*rep): 

                ind_reset = int(time_sample(1./T[t_ind[t]])/t_l_sample)

                # If the drawn time for host death is larger than the simulated time skipped it

                if ind_reset > 1E6: print('skipped: ', ind_reset)

                # Store the drawn time
                
                else: 

                    host_db_time.append(ind_reset)

                    i += 1
            
            # Create a list of hosts to die according to those times, randomly
            
            host_id = np.zeros(int(n_hosts*rep), dtype = np.int)
            
            for r in range(int(rep)):

                seq = np.arange(n_hosts)

                np.random.shuffle(seq)

                host_id[int(n_hosts*r):int(n_hosts*(1+r))] = seq
            
            # Start array to store subsamples of the microbiome of a host population
            
            subsample = np.zeros((int(n_microbe_sp+1), int(n_hosts*rep)))
            
            # Load specific simulation data
            
            data = np.load(data_dir+'data.npz')

            taxa = data['taxa']

            empty = data['empty']
            
            # Subsamples of the microbiome of a host population
            
            subsample[0, :] = empty[host_db_time, host_id]

            for j in range(1, int(n_microbe_sp+1)): subsample[j,:] = taxa[j-1][host_db_time, host_id]
            
            # Load the corresponding numeric solution
            
            p_num[t] = num_p[t_ind[t], m_ind[m]]
            
            # Compute the relevant probability from the simulations

            if c_p == 'p0': p_sim[t] = sum(subsample[ind,:] == 0)/(n_hosts*rep)

            if c_p == 'p1N': p_sim[t] = sum(subsample[ind,:] != 0)/(n_hosts*rep)

            if c_p == 'pN': p_sim[t] = sum(subsample[ind,:] == n_microbes)/(n_hosts*rep)

            print(r'$\tau$: %1.0e, $m$: %1.0e'%(T[t_ind[t]], M[m_ind[m]]))
        
        # Plot simulation data
        
        ax.plot(T[t_ind], p_sim, '^', color = colors[m])

        # Plot numeric-solution data

        ax.plot(T[t_ind], p_num, color = colors[m], linestyle = '-', label = r'%1.0e'%M[m_ind[m]])
        
        # Collect garbage
        
        gc.collect()
    
    # Set x axis log10 scale

    ax.set_xscale('log')
    
    # Set x axis limits

    ax.set_ylim(-0.02, 1.02)

    # Set title, panel label and axes labels

    ax.set_title(title_text, fontsize = 16, fontweight = 'bold', pad = 10)

    ax.text(-0.05, 1.05, label_text, fontsize = 25, fontweight = 'bold', transform = ax.transAxes)

    if y_axis == True: ax.set_ylabel(y_label, fontsize = 16)

    if y_axis != True: ax.yaxis.set_ticklabels([])
    
    ax.set_xlabel(r'prob. of host death ($\tau$)', fontsize = 16)
    
    return ax


# Figure of the probability of colonization vs. the abundance in the pool of colonizers (p)

def fig_pvspcol(ax, data_dir, title_text, label_text, y_axis, y_label):
    
    # Load data

    num_p = np.loadtxt('%sp1N.txt'%data_dir, delimiter=',')
    
    P_min = np.loadtxt('%sPmin.txt'%data_dir)

    P_max = np.loadtxt('%sPmax.txt'%data_dir)

    # Define labels and colors

    labels = [r'$10^{-7}$', r'$10^{-6}$', r'$10^{-5}$', r'$10^{-4}$', r'$10^{-3}$']

    colors = ['brown', 'red', 'gold', 'green', 'blue']
    
    # Create ranges of the frequency in the pool of colonizers 'p'
    
    interval_p = 10

    P = []

    for i in np.arange(P_min, P_max + 1): P.extend(np.arange(1, 10, 10 / interval_p) * np.power(10., i))
    
    # Plot the lines according to the values of host death rate 't'
    
    for i in range(np.shape(num_p)[0]): ax.semilogx(P, num_p[i,:], color = colors[i], label = labels[i])
    
    # Plot line indicating the maximum probability, 1
    
    ax.axhline(1, linestyle = '--', color = 'black')
    
    # Set ticks and tick labels
    
    ax.set_xticks([1E-5, 1E-4, 1E-3, 1E-2, 1E-1])

    ax.xaxis.set_ticklabels([r'$10^{-5}$', r'$10^{-4}$', r'$10^{-3}$', r'$10^{-2}$', r'$10^{-1}$'])

    # Set axes limits

    ax.set_xlim(1E-5, 1)

    ax.set_ylim(-0.02, 1.02)

    # Set title, panel label and axes labels

    ax.set_title(title_text, fontsize = 16, fontweight = 'bold', pad = 10)

    ax.text(-0.05, 1.05, label_text, fontsize = 25, fontweight = 'bold', transform = ax.transAxes)

    if y_axis == True: ax.set_ylabel(y_label, fontsize = 16)

    if y_axis != True: ax.yaxis.set_ticklabels([])
    
    ax.set_xlabel(r'freq. in the pool of colonizers ($p_1$)', fontsize = 16)
    
    return ax


# Figure of the transition between distributions shapes

def fig_dist(ax, title_text, label_text, y_axis, data):

    # Load data
    
    N = np.loadtxt(data+'N.txt')

    M_min, M_max = np.loadtxt(data+'Mmin.txt'), np.loadtxt(data+'Mmax.txt')

    T_min, T_max = np.loadtxt(data+'Tmin.txt'), np.loadtxt(data+'Tmax.txt')

    # Create ranges of microbial frequency (x_i) and host death 't' rate

    M, T = [], []

    for i in np.arange(M_min, M_max+1): M.extend(np.arange(1, 10, 2)*np.power(10., i))

    for i in np.arange(T_min, T_max+1): T.extend(np.arange(1, 10, 2)*np.power(10., i))
    
    m_ind = 20 # equivalent to m = 1E-3
    
    n = np.arange(N + 1) / N
    
    # Load probability densities
    
    Z_Phi = np.zeros((len(T), int(N + 1)))
    
    for t_ind in range(len(T)):
        
        Z_Phi[t_ind, :] = np.loadtxt(data+'/Phi/%i_%i.txt'%(m_ind, t_ind))
    
    # Set the values of probability smaller than a threshold, equal to zero

    Z_Phi[Z_Phi < 1E-9] = 1E-9
    
    # Compute the log10 of the data
    
    Z_Phi = np.log10(np.abs(Z_Phi))     

    # Create 2D space

    n_, T_ = np.meshgrid(n[1:], T)       

    # Set log10 scale of axes

    ax.set_xscale('log')

    ax.set_yscale('log')
    
    # Plot the difference btw. distributions

    ax.pcolormesh(n_, T_, Z_Phi[:,1:], cmap = mp.cm.jet, linewidth = 1E-8, zorder = 0, vmin = -9, vmax = 0)

    # Plot lines to indicate transitions of modality
    
    modetype = np.loadtxt(data+'modetype.txt', delimiter=',')[:, m_ind]
    
    for mode in np.unique(modetype):
                
        ax.axhline(T[np.where(modetype == mode)[0][-1]], color = 'k', linewidth = 2.)
    
    # Set the x tick labels

    ax.set_xticklabels(['', r'$10^{-4}$',r'$10^{-3}$',r'$10^{-2}$',r'$10^{-1}$',''])

    # Set title, panel label and axes labels

    ax.set_title(title_text, fontsize = 16, fontweight = 'bold', pad = 10)

    ax.text(-0.03, 1.03, label_text, fontsize = 25, fontweight = 'bold', transform = ax.transAxes)

    if y_axis == True: ax.set_ylabel(r'prob. of host death ($\tau$)', fontsize = 16)

    if y_axis != True: ax.yaxis.set_ticklabels([])

    ax.set_xlabel(r'freq. in hosts ($x_1$)', fontsize = 16)
        
    return ax
