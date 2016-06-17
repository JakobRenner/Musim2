#!/usr/bin/env python
# -*- coding: utf-8 -*-
#

import os
import numpy as np

import matplotlib.pyplot as plt
from matplotlib import cm

import oofem_utils as oofem

def main():

	lx,ly = 3., 5.
	nodes_x,nodes_y = 10,10
	dv = ly/2.  # Vorgeschriebene Verschiebungen oben

	infile  = 'task4.in'
	outfile = 'task4.out'

	# ==================================================================
	# ======= TODO: Inputfile-Function =================================
	# ==================================================================
	# Copy your script from task#2 for writing the inputfile into
	# a function in the file oofem_utils.py.
	# The function should have the name "write_LinElast_inputfile"
	# and should take the arguments lx,ly,nodes_x,nodes_y,dv, infile, outfile
	#
	# - 'outfile' should now show up in the first line of "...append(..) !" (I did this already...finished)
	# - dv is the vertical displacement of the top nodes
	oofem.write_LinElast_inputfile(lx,ly,nodes_x,nodes_y,dv, infile, outfile)
	# ==================================================================



	# ------- Ausfuehren von oofem aus Python heraus -------------------
	os.system("oofem -f "+infile)
	# ------------------------------------------------------------------



	# ------- Auslesen der Ergebnisse aus dem Outputfile ---------------
	# Die Daten befinden sich jetzt in 'sim'. Informationen zu der Daten-
	# struktur finden Sie in der Funktion in oofem_utils.py
	sim = oofem.readOutputfile(outfile)
	# ------------------------------------------------------------------


	# ------- test --------------------------------------------
	n=2
	ele =2
	print("average strain eps_xx at step",n,"element",ele,"=", sim.step[n].element[ele].avstrain.xx)
	# ---------------------------------------------------------



	# ==================================================================
	# ======= TODO: postprocessing and visualization ===================
	# ==================================================================
	# 1. create a numpy array 'eps' of dimension nodes_x, nodes_y
	# 2. copy the average strain avstrain.yy values from 'sim' for the last load step
	#    into the 2D-array eps
	# 3. plot the 2D array with matplotlib and compare it to the
	#    ParaView output (the .vtu file from oofem)
	#    (hint: for testing purposes it might be helpful to deform one corner
	#    of the system more - e.g. add a horizontal displacement BC to the
	#    upper left corner, compare the avstrain.yy with the paraview IST_StrainTensor 4.
	#    Thereby, you can test if your coordinate system
	#    is OK and not rotated! You can do this *directily* in the oofem input file.
	#    You should comment then the line oofem.write_LinElast_...)
	#
	# --> The next tutorial will use these data for more analysis

	# ==========================================================
	# ==================================================================
	# ======= TODO: postprocessing and visualization ===================
	# ==================================================================
	# 1. create a numpy array 'eps' of dimension nodes_x, nodes_y
	# 2. copy the average strain avstrain.yy values from 'sim' for the last load step
	#    into the 2D-array eps
	# 3. plot the 2D array with matplotlib and compare it to the
	#    ParaView output (the .vtu file from oofem)
	#    (hint: for testing purposes it might be helpful to deform one corner
	#    of the system more - e.g. add a horizontal displacement BC to the
	#    upper left corner, compare the avstrain.yy with the paraview IST_StrainTensor 4.
	#    Thereby, you can test if your coordinate system
	#    is OK and not rotated! You can do this *directily* in the oofem input file.
	#    You should comment then the line oofem.write_LinElast_...)
	#
	# --> The next tutorial will use these data for more analysis

	eps = np.zeros((nodes_y-1,nodes_x-1),float)
	X   = np.zeros((nodes_y-1,nodes_x-1),float)
	Y   = np.zeros((nodes_y-1,nodes_x-1),float)
	dx,dy = lx/(nodes_x-1.), ly/(nodes_y-1.)

	for ney in range(nodes_y-1):
		for nex in range(nodes_x-1):
			n=nex+ney*(nodes_x-1)
			eps[ney,nex] = sim.step[5].element[n].avstrain.yy
			Y[ney,nex] = dy/2. + dy*(ney)
			X[ney,nex] = dx/2. + dx*(nex)


	fig = plt.figure(figsize=(5,7),facecolor='white')
	ax1 = fig.add_subplot(111)
	plt.subplots_adjust(left=0.1, right=0.85, top=0.86, bottom=0.1, hspace = 0.3, wspace = 0.4)
	#cf1 = ax1.contourf(X, Y, eps, cmap=cm.jet)#, norm=norm, origin='lower')
	#cf1 = ax1.imshow(eps, cmap=plt.cm.jet, interpolation='nearest')
	cf1 = ax1.pcolormesh(X, Y, eps, cmap=cm.jet)#, norm=norm, origin='lower')
	fig.colorbar( cf1, orientation='vertical',aspect=20,format='%1.2e')
	ax1.grid(True)
	plt.draw()
	plt.show()
		# ==========================================================

	return 0




if __name__ == '__main__':
	main()
