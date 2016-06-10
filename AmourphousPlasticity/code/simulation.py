# -*- coding: utf-8 -*-
"""
Created on

@author: kubusRZ
"""

#! /usr/bin/env python
import greens_function
import kmc

import numpy as np

import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.colors import LogNorm
from matplotlib import ticker

def plot_KMC(time,strain,strain_pattern):

    fig, axs = plt.subplots(nrows=1, ncols=2,figsize=(12, 6), facecolor='white')
    plt.tight_layout(pad=2., w_pad=0.05)

    axs[0].plot(time,strain,'-', linewidth=2)
    axs[0].grid(True)
    axs[0].set_xlabel('Time')
    axs[0].set_ylabel('Plastic strain')

    im = axs[1].imshow(strain_pattern,cmap='afmhot')
    axs[1].set_xlabel('Plastic strain pattern')
    divider = make_axes_locatable(axs[1])
    cbar_ax1 = divider.append_axes("right", size="5%", pad=0.1)
    fig.colorbar(im, cax=cbar_ax1)

    plt.show()

L = 64
C = 0.05
max_plastic_strain = 3.0
external_stress = 0.1
P = 30.

# calculate Greens function and store it
G = greens_function.get_periodic_kernel(L,2)

# Allocate memory vectors and data to represent the system
N = L * L
stress = np.zeros(N)
plastic_strain = np.zeros(N)
energy_barrier = np.random.uniform(0., 1.0, N)
time = 0.
av_plastic_strain = 0.0

attempt_frequency = 1e12 * np.exp(-P)
#create som list to sa ve the datea that we want to plot at the end
time_data = list()
av_plastic_strain_data = list()


#test everything so far
#element = 500 # this is a munual choosen element
# now we let KMC decice which element we choose
stress += external_stress

# for _ in range(1000)
while max_plastic_strain > av_plastic_strain:
    # calculate the list of accumulated probabilties
    rates = np.exp( -P * (energy_barrier - stress) )
    probability = rates/sum(rates)
    accumulated_probabilty = np.cumsum(probability)

    r = np.random.uniform(0., 1.0)
    element = 0
    while accumulated_probabilty[element] < r:
        element += 1

    time += -np.log(r)/sum(rates)

    plastic_strain[element] += 1
    stress += C * greens_function.get_inclusion_stress_1D(G, element)
    energy_barrier[element] = np.random.uniform(0.0, 1.0) * np.exp( -plastic_strain[element])

    av_plastic_strain = np.average(plastic_strain)
    print(time, av_plastic_strain)
    time_data.append(time)
    av_plastic_strain_data.append(av_plastic_strain)

time = np.array(time)
time /= attempt_frequency
# the plotting, first with the matplotlib

#plt.imshow(stress.reshape(L, L), cmap='afmhot')
#plt.imshow(plastic_strain.reshape(L, L), cmap='afmhot')
#plt.colorbar()
#plt.show()

# fig1 = plt.figure(figsize =(7,7), facecolor = 'white')
# ax1 = fig1.add_subplot(111)
# im = ax1.imshow(abs(stress.reshape(L,L)),
#                 norm=LogNorm(vmin=1e-7, vmax = stress.max()), cmap = 'afmhot' )
# #im = ax1.imshow(plastic_strain.reshape(L,L),
# #                norm=LogNorm(vmin=1e-7, vmax = plastic_strain.max()), cmap = 'afmhot' )
# plt.colorbar(im)
# plt.show()\

# and now with the plot function
plot_KMC(time_data, av_plastic_strain_data, plastic_strain.reshape(L,L))
