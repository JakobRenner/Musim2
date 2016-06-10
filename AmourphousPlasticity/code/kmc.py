# -*- coding: utf-8 -*-
"""
Created on Thu Jan 21 14:28:54 2016

@author: dfernandez
"""

#! /usr/bin/env python

import numpy as np
import greens_function
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.colors import LogNorm
from matplotlib import ticker
import matplotlib.pyplot as plt
import numpy as np


def simulation_KMC(parameters):
    
    L = parameters[0]
    C = parameters[1]
    max_plastic_strain = parameters[2]
    external_stress = parameters[3]
    P = parameters[4]
    
    print 'Calculating Green\'s function...'
    G = greens_function.get_periodic_kernel(L,2)    
    print 'done'

    #---- allocate vectors in memory ----    
    N = L*L
    stress = np.zeros(N)
    plastic_strain = np.zeros(N)
    energy_barrier  = np.random.uniform(0.,1.0,N)
    time = 0.
    av_plastic_strain = 0.

    
    time_data = []
    av_plastic_strain_data = []
    
    
    #---- apply a load on the system ----
    stress += external_stress 

    while av_plastic_strain < max_plastic_strain: #time < max_time:#
    
        
        #---- KMC choose event ----
        rates = np.exp(-P*(energy_barrier-stress))
        probability = rates/sum(rates)
        acumulated_probabilitiy = np.cumsum(probability)
        
        r = np.random.uniform(0.,1.)
        
        element = 0
        while acumulated_probabilitiy[element] < r:
            element += 1
        
        #----- make the event happen ----
        time += -np.log(r)/sum(rates)
        plastic_strain[element] += 1
        energy_barrier[element] = np.random.uniform(0.,1.0)#*np.exp(-plastic_strain[element]/f)

        stress += C*greens_function.get_inclusion_stress_1D(G,element)
        
        #---- output data --------
        av_plastic_strain = np.average(plastic_strain)
        time_data.append(time)
        av_plastic_strain_data.append(av_plastic_strain)
        
        print time, av_plastic_strain, element

    plot_KMC(time_data,av_plastic_strain_data,plastic_strain.reshape(L,L))
    
    

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



def main():
    
    np.random.seed(128456) 

    L = 32
    C = 0.05
    max_plastic_strain = 10.
    external_stress = 0.5
    P = 50.
    parameters = [L,C,max_plastic_strain,external_stress,P]
    simulation_KMC(parameters)
    
    return 0




if __name__ == '__main__':
    main()