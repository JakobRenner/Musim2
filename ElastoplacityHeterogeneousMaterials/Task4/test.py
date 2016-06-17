#!/usr/bin/env python
# -*- coding: utf-8 -*-
#

import os
import numpy as np
import subprocess
import matplotlib.pyplot as plt
from matplotlib import cm

import oofem_utils as oofem




def main():

    lx = 3.0; ly = 4.0
    nodes_x = 10; nodes_y = 10
    dv = 0.01

    infile_elast = 'linear_elastic.in'
    outfile_elast = 'linear_elastic.out'

    infile_j2plast = 'j2_plasticity.in'
    outfile_j2plast = 'j2_plasticity.out'

    # testing the inputfiles
    oofem.write_LinElast_inputfile(lx,ly,nodes_x,nodes_y, dv, infile_elast, outfile_elast)
    oofem.write_J2plastic_inputfile(lx,ly,nodes_x,nodes_y, dv, infile_j2plast, outfile_j2plast)

    # make a directory for the output
    #os.makedirs('Output', exist_ok=True)
    #subprocess.call(["cd", "Output/"])

    # starting the simulations
    subprocess.call(["oofem", "-f", infile_elast])
    data_elast = oofem.readOutputfile(outfile_elast,True)

    avstress_yy_elast = np.zeros(len(data_elast.step))
    avstrain_yy_elast = avstress_yy_elast

    for i in range(len(data_elast.step)):
        avstress_yy_elast[i] = data_elast.step[i].element[1].avstress.yy
        avstrain_yy_elast[i] = data_elast.step[i].element[1].avstrain.yy


    subprocess.call(["oofem", "-f", infile_j2plast])

    '''
    data_j2plast = oofem.readOutputfile(outfile_j2plast)
    avstress_yy_j2plast = np.zeros(len(data_j2plast.step))
    avstrain_yy_j2plast = avstress_yy_j2plast


    for i in range(len(data_elast.step)):
        avstress_yy_j2plast[i] = data_elast.step[i].element[1].avstress.yy
        avstrain_yy_j2plast[i] = data_elast.step[i].element[1].avstrain.yy
    '''

    # the plotting
    fig, ax = plt.subplots()
    ax.plot(avstrain_yy_elast, avstress_yy_elast, marker='o', color='green', label="linear elastic")
    #ax.plot(avstrain_yy_j2plast, avstress_yy_j2plast, marker='o', color='green', label="j2 plasticity")

    plt.xlabel("Average strain_yy", size=20)
    plt.ylabel("Average stress_yy", size=20)
    #plt.xlim(-0.07, 0.07)

    ax.legend(frameon=False)
    plt.tight_layout(pad=0.5)

    #Name = ''
    #ImageFormat = 'pdf'

    #plt.savefig(Name + ImageFormat, format=ImageFormat)
    plt.show()


if __name__ == "__main__":

    main()
