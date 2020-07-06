#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@project: Stochastic colonization of microbe-free hosts (execute - figures)
@author: Román Zapién-Campos
(MPI for Evolutionary Biology - zapien@evolbio.mpg.de)
"""


# Import the relevant packages

from sc_figure import *

from matplotlib import colorbar, colors

from matplotlib.pyplot import cm



# Figure 2 (Simulations of colonization)

# Initialize figure

fig, axes = mp.subplots(nrows = 1, ncols = 2, figsize = (8.5, 5))

# Compute figure

fig_simulation(axes[0],'../data/simulations/m1.0e-02_a0.0e+00/', r'slow colonization ($\mathbf{\alpha_0 = 0}$)', 'A', True, 1E-5, r'$10^{-5}$')

fig_simulation(axes[1],'../data/simulations/m1.0e-02_a-1.0e+00/', r'fast colonization ($\mathbf{\alpha_0 \to -1}$)', 'B', False, 1E-5, r'$10^{-5}$')

# Set aspect ratio of axes

axes[0].set_aspect(1)

axes[1].set_aspect(1)

# Include colorbar to indicate the timescale

fig.subplots_adjust(bottom = 0.17, top = 1.0, left = 0.0, right = 1.0, wspace = 0.03, hspace = 0.03)

cb_ax = fig.add_axes([0.0, 0.0, 1.0, 0.04])

cb_ax.axis('off')

img = cb_ax.imshow(np.array([[0,1]]), cmap = mp.cm.gnuplot)

img.set_visible(False)

cb_ax.set_aspect('auto')

cbar = fig.colorbar(img, orientation = "horizontal", ax = cb_ax, fraction = 1.0)

cbar.ax.tick_params(labelsize = 12)

cbar.ax.set_xticklabels(cbar.ax.get_xticklabels(), rotation = 90)

cbar.set_label('Simulation timesteps ($\cdot 10^7$)',fontsize = 16)

# Save figure

mp.savefig('fig2.pdf', dpi = 300, bbox_inches = 'tight', format = 'pdf')



# Figure 3 (Microbial load)

# Indicate contours

c_levels = [-4, -3, -2, - 1, -5.23E-1, -1E-1, -1E-2, -1E-5]

c_levels_labels = {-4:r'$10^{-4}$', -3:r'$10^{-3}$', -2:r'$0.01$', -1:r'$0.1$', -5.23E-1:r'0.3', -1E-1:r'0.8', -1E-2:r'$0.98$', -1E-5:r'$\approx 1$'}

# Define colors for the 'modality' section of the figure 

cmap = cm.Set1

cmaplist_plot = [cmap(1), cmap(0), cmap(2), cmap(3), cmap(4)]

cmap_plot = colors.LinearSegmentedColormap.from_list('custom cmap', cmaplist_plot, 5)

cmaplist_legend = [cmap(0), cmap(1), cmap(2), cmap(3), cmap(4)]

cmap_legend = colors.LinearSegmentedColormap.from_list('custom cmap', cmaplist_legend, 5)

bounds = np.linspace(1, 5, 5)

norm = colors.BoundaryNorm(bounds, 5)

# Initialize figure

fig, axes = mp.subplots(nrows = 2, ncols = 2, figsize = (6.7, 8.0))

# Compute figure

axes[0,0] = fig_dp(axes[0,0], r'slow col. ($\mathbf{\alpha_{0} = 0.0}$)', 'A', True, '../data/numerics/N_10000_p_0._q_1._a_0./', c_levels, c_levels_labels)

axes[0,1] = fig_dp(axes[0,1], r'fast col. ($\mathbf{\alpha_{0} = -0.9}$)', 'B', False, '../data/numerics/N_10000_p_0._q_1._a_-0.9/', c_levels, c_levels_labels)

axes[1,0] = fig_mode(axes[1,0], r'slow col. ($\mathbf{\alpha_{0} = 0.0}$)', 'C', True, '../data/numerics/N_10000_p_0._q_1._a_0./', cmap_plot, norm)

axes[1,1] = fig_mode(axes[1,1], r'fast col. ($\mathbf{\alpha_{0} = -0.9}$)', 'D', False, '../data/numerics/N_10000_p_0._q_1._a_-0.9/', cmap_plot, norm)

# Set aspect ratio of axes

axes[0,0].set_aspect(35./40)

axes[0,1].set_aspect(35./40)

axes[1,0].set_aspect(35./40)

axes[1,1].set_aspect(35./40)

# Set separation between panels

fig.subplots_adjust(bottom = 0.0, top = 1.0, left = 0.0, right = 1.0, wspace = 0.03, hspace = 0.44)

# Add secondary titles

mp.text(0.5, 1.07, 'Difference between distributions', fontsize = 16, transform = fig.transFigure, horizontalalignment = 'center', weight = 'bold')

mp.text(0.5, 0.48, 'Shape of stationary distribution', fontsize = 16, transform = fig.transFigure, horizontalalignment = 'center', weight = 'bold')

# Save figure

mp.savefig('fig3.pdf', dpi = 300, bbox_inches = 'tight', format = 'pdf')



# Figure 4 (Microbial taxon 1)

# Indicate contours

c_levels = [-4, -3, -2, -1, -5.23E-1, -1E-1, -1E-2, -1E-5]

c_levels_labels = {-4:r'$10^{-4}$', -3:r'$10^{-3}$', -2:r'$0.01$', -1:r'$0.1$', -5.23E-1:r'0.3', -1E-1:r'0.8', -1E-2:r'$0.98$', -1E-5:r'$\approx 1$'}

# Define colors for the 'modality' section of the figure

cmap = cm.Set1

cmaplist_plot = [cmap(0), cmap(1), cmap(2), cmap(3), cmap(4)]

cmap_plot = colors.LinearSegmentedColormap.from_list('custom cmap', cmaplist_plot, 5)

cmaplist_legend = [cmap(0), cmap(1), cmap(2), cmap(3), cmap(4)]

cmap_legend = colors.LinearSegmentedColormap.from_list('custom cmap', cmaplist_legend, 5)

bounds = np.linspace(1, 5, 5)

norm = colors.BoundaryNorm(bounds, 5)

# Initialize figure

fig, axes = mp.subplots(nrows = 2, ncols = 3, figsize = (8.5, 8.0))

# Compute figure

axes[0,0] = fig_dp(axes[0,0], r'$\mathbf{p_1 = 1.0}$', 'A', True, '../data/numerics/N_10000_p_1._q_0._a_0./', c_levels, c_levels_labels)

axes[0,1] = fig_dp(axes[0,1], r'$\mathbf{p_1 = 0.5}$', 'B', False, '../data/numerics/N_10000_p_0.5_q_0._a_0./', c_levels, c_levels_labels)

axes[0,2] = fig_dp(axes[0,2], r'$\mathbf{p_1 = 0.1}$', 'C', False, '../data/numerics/N_10000_p_0.1_q_0._a_0./', c_levels, c_levels_labels)

axes[1,0] = fig_mode(axes[1,0], r'$\mathbf{p_1 = 1.0}$', 'D', True, '../data/numerics/N_10000_p_1._q_0._a_0./', cmap_plot, norm)

axes[1,1] = fig_mode(axes[1,1], r'$\mathbf{p_1 = 0.5}$', 'E', False, '../data/numerics/N_10000_p_0.5_q_0._a_0./', cmap_plot, norm)

axes[1,2] = fig_mode(axes[1,2], r'$\mathbf{p_1 = 0.1}$', 'F', False, '../data/numerics/N_10000_p_0.1_q_0._a_0./', cmap_plot, norm)

# Set aspect ratio of axes

axes[0,0].set_aspect(35./40)

axes[0,1].set_aspect(35./40)

axes[0,2].set_aspect(35./40)

axes[1,0].set_aspect(35./40)

axes[1,1].set_aspect(35./40)

axes[1,2].set_aspect(35./40)

# Set separation between panels

fig.subplots_adjust(bottom = 0.0, top = 1.0, left = 0.0, right = 1.0, wspace = 0.03, hspace = 0.1)

# Add secondary titles

mp.text(0.5, 1.0, 'Difference between distributions', fontsize = 16, transform = fig.transFigure, horizontalalignment = 'center', weight = 'bold')

mp.text(0.5, 0.48, 'Shape of stationary distribution', fontsize = 16, transform = fig.transFigure, horizontalalignment = 'center', weight = 'bold')

# Save figure

mp.savefig('fig4.pdf', dpi = 300, bbox_inches = 'tight', format = 'pdf')



# Figure 5 (Probability of full colonization)

# Indicate the observable

c_p = 'p0'

# Indicate contours

c_levels = [-5, -2, -1, -3E-1, -1E-1, -1E-2, -1E-3, -1E-5]

c_levels_labels = {-5:r'$10^{-5}$', -2:'0.01', -1:'0.1', -3E-1:r'0.5', -1E-1:r'0.8', -1E-2:r'$0.98$', -1E-3:r'$0.998$', -1E-5:r'$\approx 1$'}

# Initialize figure

fig, axes = mp.subplots(nrows = 1, ncols = 2, figsize = (8, 4))

# Compute figure

axes[0] = fig_p(axes[0], r'slow col. ($\mathbf{\alpha_{0} = 0.0}$)', 'A', True, '../data/numerics/N_10000_p_0._q_1._a_0./', c_p, c_levels, c_levels_labels)

axes[1] = fig_p(axes[1], r'fast col. ($\mathbf{\alpha_{0} = -0.9}$)', 'B', False, '../data/numerics/N_10000_p_0._q_1._a_-0.9/', c_p, c_levels, c_levels_labels)

# Set aspect ratio of axes

axes[0].set_aspect(35./40)

axes[1].set_aspect(35./40)

# Set separation between panels

fig.subplots_adjust(bottom = 0.0, top = 1.0, left = 0.0, right = 1.0, wspace = 0.03, hspace = 0.02)

# Save figure

mp.savefig('fig5.pdf', dpi = 300, bbox_inches = 'tight', format = 'pdf')



# Figure 6 (Probability of colonization of taxon 1)

# Indicate the observable

c_p = 'p1N'

# Indicate contours

c_levels = [-4, -2, -1, -5.23E-1, -3E-1, -1E-1, -1E-2, -1E-3, -1E-5]

c_levels_labels = {-4:r'$10^{-4}$', -2:r'0.01', -1:r'0.1', -5.23E-1:r'0.3', -3E-1:r'0.5', -1E-1:r'0.8', -1E-2:r'$0.98$', -1E-3:r'$0.998$', -1E-5:r'$\approx 1$'}

# Initialize figure

fig, axes = mp.subplots(nrows = 1, ncols = 3, figsize = (8.5, 5))

# Compute figure

axes[0] = fig_p(axes[0], r'$\mathbf{p_1 = 1.0}$', 'A', True, '../data/numerics/N_10000_p_1._q_0._a_0./', c_p, c_levels, c_levels_labels)

axes[1] = fig_p(axes[1], r'$\mathbf{p_1 = 0.5}$', 'B', False, '../data/numerics/N_10000_p_0.5_q_0._a_0./', c_p, c_levels, c_levels_labels)

axes[2] = fig_p(axes[2], r'$\mathbf{p_1 = 0.1}$', 'C', False, '../data/numerics/N_10000_p_0.1_q_0._a_0./', c_p, c_levels, c_levels_labels)

# Set aspect ratio of axes

axes[0].set_aspect(35./40)

axes[1].set_aspect(35./40)

axes[2].set_aspect(35./40)

# Set separation between panels

fig.subplots_adjust(bottom = 0.0, top = 1.0, left = 0.0, right = 1.0, wspace = 0.03, hspace = 0.02)

# Save figure

mp.savefig('fig6.pdf', dpi = 300, bbox_inches = 'tight', format = 'pdf')



# Figure Sup. 1 (Simulations of colonization: limited migration)

# Initialize figure

fig, axes = mp.subplots(nrows = 1, ncols = 2, figsize = (8.5, 5))

# Compute figure

fig_simulation(axes[0],'../data/simulations/m1.0e-04_a0.0e+00/', r'slow colonization ($\mathbf{\alpha_0 = 0})$', 'A', True, 1E-5, r'$10^{-5}$')

fig_simulation(axes[1],'../data/simulations/m1.0e-04_a-1.0e+00/', r'fast colonization ($\mathbf{\alpha_0 \to -1}$)', 'B', False, 1E-5, r'$10^{-5}$')

# Set aspect ratio of axes

axes[0].set_aspect(1)

axes[1].set_aspect(1)

# Include colorbar to indicate the timescale

fig.subplots_adjust(bottom = 0.17, top = 1.0, left = 0.0, right = 1.0, wspace = 0.03, hspace = 0.03)

cb_ax = fig.add_axes([0.0, 0.0, 1.0, 0.04])

cb_ax.axis('off')

img = cb_ax.imshow(np.array([[0,1]]), cmap = mp.cm.gnuplot)

img.set_visible(False)

cb_ax.set_aspect('auto')

cbar = fig.colorbar(img, orientation = "horizontal", ax = cb_ax, fraction = 1.0)

cbar.ax.tick_params(labelsize=12)

cbar.ax.set_xticklabels(cbar.ax.get_xticklabels(), rotation = 90)

cbar.set_label('Simulation timesteps ($\cdot 10^7$)',fontsize = 16)

# Save figure

mp.savefig('figsup1.pdf', dpi = 300, bbox_inches = 'tight', format = 'pdf')



# Figure Sup. 2 (Probability of full colonization: smaller microbiome size)

# Indicate the observable

c_p = 'p0'

# Indicate contours

c_levels = [-5, -2, -1, -3E-1, -1E-1, -1E-2, -1E-3, -1E-5]

c_levels_labels = {-5:r'$10^{-5}$', -2:'0.01', -1:'0.1', -3E-1:r'0.5', -1E-1:r'0.8', -1E-2:r'$0.98$', -1E-3:r'$0.998$', -1E-5:r'$\approx 1$'}

# Initialize figure

fig, axes = mp.subplots(nrows = 1, ncols = 2, figsize = (8, 4))

# Compute figure

axes[0] = fig_p(axes[0], r'slow col. ($\mathbf{\alpha_{0} = 0.0}$)', 'A', True, '../data/numerics/N_1000_p_0._q_1._a_0./', c_p, c_levels, c_levels_labels)

axes[1] = fig_p(axes[1], r'fast col. ($\mathbf{\alpha_{0} = -0.9}$)', 'B', False, '../data/numerics/N_1000_p_0._q_1._a_-0.9/', c_p, c_levels, c_levels_labels)

# Set aspect ratio of axes

axes[0].set_aspect(35./40)

axes[1].set_aspect(35./40)

# Set separation between panels

fig.subplots_adjust(bottom = 0.0, top = 1.0, left = 0.0, right = 1.0, wspace = 0.03, hspace = 0.02)

# Save figure

mp.savefig('figsup2.pdf', dpi = 300, bbox_inches = 'tight', format = 'pdf')



# Figure Sup. 3 (Probability of microbe-free hosts)

# Indicate the observable

c_p = 'pN'

# Indicate contours

c_levels = [-5, -3, -2, -1, -3E-1, -1E-1, -1E-2, -1E-3]

c_levels_labels = {-5:r'$10^{-5}$', -3:r'$10^{-3}$', -2:r'0.01', -1:r'0.1', -5.23E-1:r'0.3', -3E-1:r'0.5', -1E-1:r'0.8', -1E-2:r'$0.98$', -1E-3:r'$0.998$'}

# Initialize figure

fig, axes = mp.subplots(nrows = 1, ncols = 2, figsize = (8, 4))

# Compute figure

axes[0] = fig_p(axes[0], r'slow col. ($\mathbf{\alpha_{0} = 0.0}$)', 'A', True, '../data/numerics/N_10000_p_0._q_1._a_0./', c_p, c_levels, c_levels_labels)

axes[1] = fig_p(axes[1], r'fast col. ($\mathbf{\alpha_{0} = -0.9}$)', 'B', False, '../data/numerics/N_10000_p_0._q_1._a_-0.9/', c_p, c_levels, c_levels_labels)

# Set aspect ratio of axes

axes[0].set_aspect(35./40)

axes[1].set_aspect(35./40)

# Set separation between panels

fig.subplots_adjust(bottom = 0.0, top = 1.0, left = 0.0, right = 1.0, wspace = 0.03, hspace = 0.02)

# Save figure

mp.savefig('figsup3.pdf', dpi = 300, bbox_inches = 'tight', format = 'pdf')



# Figure Sup. 4 (Simulations vs. model: probability of microbe-free hosts)

# Indicate the observable

c_p = 'pN'

c_t = 'empty-space'

y_label = r'$P[x_0 > (N-1)/N]$'

# Initialize figure

fig, axes = mp.subplots(nrows = 1, ncols = 2, figsize = (8, 4))

# Compute figure

axes[0] = fig_simvsmodel(axes[0], r'slow  prolif. ($\mathbf{\alpha_{0} = 0.0}$)', 'A', True, '../data/numerics/N_10000_p_0._q_1._a_0./', 'a0.0e+00', c_p, c_t, y_label)

axes[1] = fig_simvsmodel(axes[1], r'fast prolif. ($\mathbf{\alpha_{0} = -0.9}$)', 'B', False, '../data/numerics/N_10000_p_0._q_1._a_-0.9/', 'a-9.0e-01', c_p, c_t, y_label)

# Add text labels including colorcode

colors = ['brown','red','orange','gold','green','blue','turquoise']

leg_labels = [r'$m = 10^{%i}$'%i for i in range(-7,0,1)]

for i in range(7): axes[0].text(1E-6, 0.2+0.05*i, leg_labels[i], ha = 'center', va = 'center', fontsize = 14, color = colors[i])

for i in range(7): axes[1].text(1E-6, 0.2+0.05*i, leg_labels[i], ha = 'center', va = 'center', fontsize = 14, color = colors[i])

# Set separation between panels

mp.subplots_adjust(bottom = 0.0, top = 1.0, left = 0.0, right = 1.0, wspace = 0.03)

# Save figure

mp.savefig('figsup4.pdf', dpi = 300, bbox_inches = 'tight', format = 'pdf')



# Figure Sup. 5 (Simulations vs. model: probability of full colonization)

# Indicate the observable

c_p = 'p0'

c_t = 'empty-space'

y_label = r'$P[x_0 < 1/N]$'

# Initialize figure

fig, axes = mp.subplots(nrows = 1, ncols = 2, figsize = (8, 4))

# Compute figure

axes[0] = fig_simvsmodel(axes[0], r'slow prolif. ($\mathbf{\alpha_{0} = 0.0}$)', 'A', True, '../data/numerics/N_10000_p_0._q_1._a_0./', 'a0.0e+00', c_p, c_t, y_label)

axes[1] = fig_simvsmodel(axes[1], r'fast prolif. ($\mathbf{\alpha_{0} = -0.9}$)', 'B', False, '../data/numerics/N_10000_p_0._q_1._a_-0.9/', 'a-9.0e-01', c_p, c_t, y_label)

# Add text labels including colorcode

colors = ['brown','red','orange','gold','green','blue','turquoise']

leg_labels = [r'$m = 10^{%i}$'%i for i in range(-7,0,1)]

for i in range(7): axes[0].text(1E-4, 0.25+0.05*i, leg_labels[i], ha = 'center', va = 'center', fontsize = 14, color = colors[i])

for i in range(7): axes[1].text(1E-4, 0.45+0.05*i, leg_labels[i], ha = 'center', va = 'center', fontsize = 14, color = colors[i])

# Set separation between panels

mp.subplots_adjust(bottom = 0.0, top = 1.0, left = 0.0, right = 1.0, wspace = 0.03)

# Save figure

mp.savefig('figsup5.pdf', dpi = 300, bbox_inches = 'tight', format = 'pdf')



# Figure Sup. 6 (Simulations vs. model: probability of colonization of taxon 1)

# Indicate the observable

c_p = 'p1N'

c_t = 'taxa'

y_label = r'$P[x_1 \geq 1/N]$'

# Initialize figure

fig, axes = mp.subplots(nrows = 1, ncols = 1, figsize = (4, 4))

# Compute figure

axes = fig_simvsmodel(axes, r'$\mathbf{p_1 = 0.5}$', '', True, '../data/numerics/N_10000_p_0.5_q_0._a_0./', 'a0.0e+00', c_p, c_t, y_label)

# Add text labels including colorcode

colors = ['brown','red','orange','gold','green','blue','turquoise']

leg_labels = [r'$m = 10^{%i}$'%i for i in range(-7,0,1)]

for i in range(7): axes.text(2E-6, 0.25+0.05*i, leg_labels[i], ha = 'center', va = 'center', fontsize = 14, color = colors[i])

# Set separation between panels

mp.subplots_adjust(bottom = 0.0, top = 1.0, left = 0.0, right = 1.0, wspace = 0.03, hspace = 0.02)

# Save figure

mp.savefig('figsup6.pdf', dpi = 300, bbox_inches = 'tight', format = 'pdf')



# Figure Sup. 7 (Probability of colonization of taxon 1 as a function of p)

# Indicate the observable

y_label = r'$P[x_1 \geq 1/N]$'

# Initialize figure

fig, axes = mp.subplots(nrows = 1, ncols = 2, figsize = (8, 4))

# Compute figure

axes[0] = fig_pvspcol(axes[0], '../data/numerics/N_10000_p_0.5_q_0._a_0._m_0.01/', r'$\mathbf{m = 0.01}$', 'A', True, y_label)

axes[1] = fig_pvspcol(axes[1], '../data/numerics/N_10000_p_0.5_q_0._a_0._m_0.1/', r'$\mathbf{m = 0.1}$', 'B', False, '')

# Set aspect ratio of axes

axes[0].set_aspect(1)

axes[1].set_aspect(1)

# Add text labels including colorcode

colors = ['brown','red','gold','green','blue']

leg_labels = [r'$\tau = 10^{%i}$'%i for i in range(-7,-2,1)]

for i in range(5): axes[0].text(5E-5, 0.4+0.1*i, leg_labels[i], ha = 'center', va = 'center', fontsize = 14, color = colors[i])

for i in range(5): axes[1].text(5E-5, 0.4+0.1*i, leg_labels[i], ha = 'center', va = 'center', fontsize = 14, color = colors[i])

# Set separation between panels

mp.subplots_adjust(bottom = 0.0, top = 1.0, left = 0.0, right = 1.0, wspace = 0.03, hspace = 0.0)

# Save figure

mp.savefig('figsup7.pdf', dpi = 300, bbox_inches = 'tight', format = 'pdf')