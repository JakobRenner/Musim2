#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  untitled.py
#
#  Copyright 2014 Stefan Sandfeld <stefan.sandfeld@faude>
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
#  MA 02110-1301, USA.
#
#

import numpy as np
import numpy.random as rnd


def main():

	# size of the system and number of nodes
	lx,ly = 3., 5.
	nodes_x,nodes_y = 30,40

	dofs = nodes_x*nodes_y
	Nele=(nodes_x-1)*(nodes_y-1)


	oofem_input = list()
	oofem_input.append('task2.in.out')
	oofem_input.append('Task2')
	oofem_input.append('NonLinearStatic nsteps 10 controllmode 1 deltaT 0.1 rtolv 1e-6  MaxIter 1000   nmodules 3')
	oofem_input.append('vtkxml tstep_all domain_all primvars 1 1 vars 4 1 2 4 5   stype 1')
	oofem_input.append('hom tstep_all')
	oofem_input.append('GPExportModule tstep_all domain_all ncoords 3 vars 3 1 4 27')
	oofem_input.append('domain PlaneStrain')
	oofem_input.append('OutputManager tstep_all dofman_all element_all')
	oofem_input.append('ndofman {} nelem {} ncrosssect 1 nmat {} nbc 2 nic 0 nltf 1'.format(dofs,Nele,Nele))

	node_id = 0

	# ------------------------------------------------------
	# ------- TODO: define nodes with 2 for-loops ----------
	# ------------------------------------------------------
	# hint: if-elif-... takes care of BCs at bottom and top
	# number of nodes can be counted by incrementing node_id
	# the 2 loops should go over all node numbers in x and y direction
	# node numbers start at 1 and not at 0 !
	# ------------------------------------------------------
	for i in range(nodes_y):
		for j in range(nodes_x):
			node_id += 1
			# the first row
			if node_id <= nodes_x:

				x = j * lx/(nodes_x-1)
				y = 0.0
				#print(node_id, x, y)
				oofem_input.append('node {} coords 3 \t {} \t {} \t 0.0   bc 2 1 1'.format(node_id, x, y))
			# the nodes inbetween
			elif node_id > nodes_x and node_id <= dofs-nodes_x:
				x = j * lx/(nodes_x-1)
				y = i * ly/(nodes_y-1)
				oofem_input.append('node {} coords 3 \t {} \t {} \t 0.0 '.format(node_id, x, y))
			# the last row
			elif node_id >= (dofs-nodes_x+1):
				x = j * lx/(nodes_x-1)
				y = i * ly/(nodes_y-1)
				oofem_input.append('node {} coords 3 \t {} \t {} \t 0.0   bc 2 0 2'.format(node_id, x, y))

	# ------------------------------------------------------
	# ------- TODO: define elements with 2 for-loops -------
	# ------------------------------------------------------
	# define: node numbers of the element (lowerleft, lower right, upper left, upper right)
	# each element should get its own material id
	# element numbers start at 1 and not at 0 !
	# ------------------------------------------------------
	element_id = 0
	for i in range(nodes_y-1):
		for j in range(nodes_x-1):
			element_id += 1
			lower_left = element_id + i
			lower_right = element_id + 1 + i
			upper_right = lower_right + nodes_x
			upper_left = lower_left + nodes_x
			oofem_input.append('quad1planestrain \t {} nodes 4 \t {} \t {} \t {} \t{} crossSect 1 mat {}'.format(element_id, lower_left, lower_right, upper_right, upper_left, element_id))
	# ------------------------------------------------------
	# ------- define CS ------------------------------------
	# ------------------------------------------------------
	oofem_input.append('SimpleCS 1 thick 1')
	# ------------------------------------------------------


	# ------------------------------------------------------
	# ------- TODO: define materials
	# ------------------------------------------------------
	# use loop over all elements
	# Young's modulus should be random
	# ------------------------------------------------------
	element_id = 0
	for i in range(nodes_y-1):
		for j in range(nodes_x-1):
			element_id += 1
			random_E = rnd.uniform(low=0.5, high=1.0)
			oofem_input.append('IsoLE {} d 0. E {} n 0.3 tAlpha 0'.format(element_id, random_E))


	# ------------------------------------------------------
	# ------- define 2 BCs ---------------------------------
	# ------------------------------------------------------
	oofem_input.append('BoundaryCondition 1 loadTimeFunction 1 prescribedvalue 0')
	oofem_input.append('BoundaryCondition 2 loadTimeFunction 1 prescribedvalue 0.002')
	# ------------------------------------------------------


	# ------------------------------------------------------
	# ------- define time dependent function ---------------
	# ------------------------------------------------------
	oofem_input.append('PiecewiseLinFunction 1 nPoints 2 t 2 0.0 1.0 f(t) 2 0.0 1.0')
	# ------------------------------------------------------



	fname = 'task2_big.in'
	with open(fname, 'w') as f:
		f.write('\n'.join(oofem_input))

	return 0



if __name__ == '__main__':
	main()
