# -*- coding: utf-8 -*-
"""
Created on Thu Jan 23 14:06:50 2014

@author: manu

"""


import numpy as np
import numpy.random as rnd
from collections import namedtuple



def write_LinElast_inputfile(lx,ly,nodes_x,nodes_y,dv, infile, outfile):
	""" lx,ly: size of the system and number of nodes
	    nodes_x,nodes_y: number of nodes into x and y direction
	    dv: prescribed displacement of the top row of nodes
	    infile, outfile: filenames for inputfile and results file
	"""
	
	dofs = nodes_x*nodes_y
	Nele=(nodes_x-1)*(nodes_y-1)
		
	oofem_input = list()
	oofem_input.append(outfile)
	oofem_input.append('Linear elastic FEM')
	oofem_input.append('NonLinearStatic nsteps 10 controllmode 1 deltaT 0.1 rtolv 1e-6  MaxIter 1000   nmodules 3')
	oofem_input.append('vtkxml tstep_all domain_all primvars 1 1 vars 4 1 2 4 5   stype 1')
	oofem_input.append('hom tstep_all')
	oofem_input.append('GPExportModule tstep_all domain_all ncoords 3 vars 3 1 4 27')
	oofem_input.append('domain PlaneStrain')
	oofem_input.append('OutputManager tstep_all dofman_all element_all')
	oofem_input.append('ndofman {} nelem {} ncrosssect 1 nmat {} nbc 2 nic 0 nltf 1'.format(dofs,Nele,Nele))

	# add your solution here! ###################################

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
			if node_id == 1:
				x = j * lx / (nodes_x - 1)
				y = 0.0
				oofem_input.append('node {} coords 3 \t {} \t {} \t 0.0   bc 2 1 1'.format(node_id, x, y))
			elif node_id <= nodes_x:
				x = j * lx / (nodes_x - 1)
				y = 0.0
				# print(node_id, x, y)
				oofem_input.append('node {} coords 3 \t {} \t {} \t 0.0   bc 2 0 1'.format(node_id, x, y))
			# the nodes inbetween
			elif node_id > nodes_x and node_id <= dofs - nodes_x:
				x = j * lx / (nodes_x - 1)
				y = i * ly / (nodes_y - 1)
				oofem_input.append('node {} coords 3 \t {} \t {} \t 0.0 '.format(node_id, x, y))
			# the last row
			elif node_id >= (dofs - nodes_x + 1):
				x = j * lx / (nodes_x - 1)
				y = i * ly / (nodes_y - 1)
				oofem_input.append('node {} coords 3 \t {} \t {} \t 0.0   bc 2 0 2'.format(node_id, x, y))

	# ------------------------------------------------------
	# ------- TODO: define elements with 2 for-loops -------
	# ------------------------------------------------------
	# define: node numbers of the element (lowerleft, lower right, upper left, upper right)
	# each element should get its own material id
	# element numbers start at 1 and not at 0 !
	# ------------------------------------------------------
	element_id = 0
	for i in range(nodes_y - 1):
		for j in range(nodes_x - 1):
			element_id += 1
			lower_left = element_id + i
			lower_right = element_id + 1 + i
			upper_right = lower_right + nodes_x
			upper_left = lower_left + nodes_x
			oofem_input.append(
				'quad1planestrain \t {} nodes 4 \t {} \t {} \t {} \t{} crossSect 1 mat {}'.format(element_id,
																								  lower_left,
																								  lower_right,
																								  upper_right,
																								  upper_left,
																								  element_id))
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
	for i in range(nodes_y - 1):
		for j in range(nodes_x - 1):
			element_id += 1
			random_E = rnd.uniform(low=0.5, high=1.0)
			oofem_input.append('IsoLE {} d 0. E {} n 0.3 tAlpha 0'.format(element_id, random_E))

	# ------------------------------------------------------
	# ------- define 2 BCs ---------------------------------
	# ------------------------------------------------------
	oofem_input.append('BoundaryCondition 1 loadTimeFunction 1 prescribedvalue 0')
	oofem_input.append('BoundaryCondition 2 loadTimeFunction 1 prescribedvalue {}'.format(dv))
	# ------------------------------------------------------


	# ------------------------------------------------------
	# ------- define time dependent function ---------------
	# ------------------------------------------------------
	oofem_input.append('PiecewiseLinFunction 1 nPoints 2 t 2 0.0 1.0 f(t) 2 0.0 1.0')
	# ------------------------------------------------------

	with open(infile, 'w') as f:
		f.write('\n'.join(oofem_input))
	
	return 0


def write_J2plasti_inputfile(lx, ly, nodes_x, nodes_y, dv, infile, outfile):
    """ lx,ly: size of the system and number of nodes
        nodes_x,nodes_y: number of nodes into x and y direction
        dv: prescribed displacement of the top row of nodes
        infile, outfile: filenames for inputfile and results file
    """

    dofs = nodes_x * nodes_y
    Nele = (nodes_x - 1) * (nodes_y - 1)

    oofem_input = list()
    oofem_input.append(outfile)
    oofem_input.append('J2 plasticity FEM')
    oofem_input.append('NonLinearStatic nsteps 10 controllmode 1 deltaT 0.1 rtolv 1e-6  MaxIter 1000   nmodules 3')
    oofem_input.append('vtkxml tstep_all domain_all primvars 1 1 vars 4 1 2 4 5   stype 1')
    oofem_input.append('hom tstep_all')
    oofem_input.append('GPExportModule tstep_all domain_all ncoords 3 vars 3 1 4 27')
    oofem_input.append('domain PlaneStrain')
    oofem_input.append('OutputManager tstep_all dofman_all element_all')
    oofem_input.append('ndofman {} nelem {} ncrosssect 1 nmat {} nbc 2 nic 0 nltf 1'.format(dofs, Nele, Nele))

    # add your solution here! ###################################

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
            if node_id == 1:
                x = j * lx / (nodes_x - 1)
                y = 0.0
                oofem_input.append('node {} coords 3 \t {} \t {} \t 0.0   bc 2 1 1'.format(node_id, x, y))
            elif node_id <= nodes_x:
                x = j * lx / (nodes_x - 1)
                y = 0.0
                # print(node_id, x, y)
                oofem_input.append('node {} coords 3 \t {} \t {} \t 0.0   bc 2 0 1'.format(node_id, x, y))
            # the nodes inbetween
            elif node_id > nodes_x and node_id <= dofs - nodes_x:
                x = j * lx / (nodes_x - 1)
                y = i * ly / (nodes_y - 1)
                oofem_input.append('node {} coords 3 \t {} \t {} \t 0.0 '.format(node_id, x, y))
            # the last row
            elif node_id >= (dofs - nodes_x + 1):
                x = j * lx / (nodes_x - 1)
                y = i * ly / (nodes_y - 1)
                oofem_input.append('node {} coords 3 \t {} \t {} \t 0.0   bc 2 0 2'.format(node_id, x, y))

    # ------------------------------------------------------
    # ------- TODO: define elements with 2 for-loops -------
    # ------------------------------------------------------
    # define: node numbers of the element (lowerleft, lower right, upper left, upper right)
    # each element should get its own material id
    # element numbers start at 1 and not at 0 !
    # ------------------------------------------------------
    element_id = 0
    for i in range(nodes_y - 1):
        for j in range(nodes_x - 1):
            element_id += 1
            lower_left = element_id + i
            lower_right = element_id + 1 + i
            upper_right = lower_right + nodes_x
            upper_left = lower_left + nodes_x
            oofem_input.append(
                'quad1planestrain \t {} nodes 4 \t {} \t {} \t {} \t{} crossSect 1 mat {}'.format(element_id,
                                                                                                  lower_left,
                                                                                                  lower_right,
                                                                                                  upper_right,
                                                                                                  upper_left,
                                                                                                  element_id))
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
    for i in range(nodes_y - 1):
        for j in range(nodes_x - 1):
            element_id += 1
            random_E = rnd.uniform(low=0.5, high=1.0)
            oofem_input.append('J2mat 2 d 1. Ry 0.3 E {} n 0.33 IHM 0.4 tAlpha 0'.format(element_id, random_E))

    # ------------------------------------------------------
    # ------- define 2 BCs ---------------------------------
    # ------------------------------------------------------
    oofem_input.append('BoundaryCondition 1 loadTimeFunction 1 prescribedvalue 0')
    oofem_input.append('BoundaryCondition 2 loadTimeFunction 1 prescribedvalue {}'.format(dv))
    # ------------------------------------------------------


    # ------------------------------------------------------
    # ------- define time dependent function ---------------
    # ------------------------------------------------------
    oofem_input.append('PiecewiseLinFunction 1 nPoints 2 t 2 0.0 1.0 f(t) 2 0.0 1.0')
    # ------------------------------------------------------

    with open(infile, 'w') as f:
        f.write('\n'.join(oofem_input))

    return 0




def readOutputfile(filename, verbose=False):
	"""ReadOut. This function reads a ooFem .out file and extracts
	the nodes displacements in 2 dimensions, and the strains and stresses 
	at 4 Gauss points of each element. This done for every step of the 
	simulation found in the file.
	The the element average of stresses and strains at Gauss points is
	calculated.
	
	The function returns a variable 'SimData' which contains all extracted 
	information. The variable has the following structure:
	
	- Get the stresses at a certain Gauss point:
		SimData.steGGqp[nStep].element[nElement].gp[nGP].stress.xx
		SimData.step[nStep].element[nElement].gp[nGP].stress.yy
		SimData.step[nStep].element[nElement].gp[nGP].stress.zz
		SimData.step[nStep].element[nElement].gp[nGP].stress.yz
		SimData.step[nStep].element[nElement].gp[nGP].stress.zx
		SimData.step[nStep].element[nElement].gp[nGP].stress.xy
		
	- Get the strains at a certain Gauss point:
		SimData.step[nStep].element[nElement].gp[nGP].strain.xx
		
	- Get the node displacement in certain dimension
		SimData.step[nStep].node[nNode].ux
		SimData.step[nStep].node[nNode].uy
		
	- Get the average stress/strains per element
		SimData.step[nStep].element[nElement].avstress.xx
		SimData.step[nStep].element[nElement].avstrain.xx
	
	
	Not implemented yet:
	
	- Get the equivalent strains per element
		SimData.step[nStep].element[nElement].eqstrain.xx

	"""	
	
	# -----------------------------------------------------------------------------
	# Defining the classes for data structure
    T_Simulation = namedtuple('Simulation', ['step'])
    T_Step = namedtuple('Step', ['element', 'node'])

    T_Displacement = namedtuple('Displacement', ['ux', 'uy'])

    T_Element = namedtuple('Element', ['gp', 'avstrain', 'avstress', 'eqstrain'])
    T_GP = namedtuple('GP', ['stress', 'strain'])
    T_Stresses = namedtuple('Stresses', ['xx', 'yy', 'zz', 'yz', 'zx', 'xy'])
    T_Strains = namedtuple('Strains', ['xx', 'yy', 'zz', 'yz', 'zx', 'xy'])
	# -----------------------------------------------------------------------------
	
	nSteps = 0 # Simulation step counter
	
	SimData = T_Simulation(list())
	
	with open(filename) as f:
		line = f.readline() # Read in the first line of the input file

		while True: # Loop over all lines of the input file

            # Read the nodes displacements
			if line == 'DofManager output:\n': # String starts a list of nodes displacement information
				nSteps += 1 # The above string starts a new simulation step
				line = f.readline() # Cancel ---------- seperator
				line = f.readline()
				Nodes = list() # Initialize/clear list of nodes

				while line != '\n' and line != 'Element output:\n': # Strings that finish the list
	#				nNode = int(line.strip().split()[1]) # Node id
					line = f.readline()
					dim1 = float(line.strip().split()[3]) # Displacement dim1
					line = f.readline()
					dim2 = float(line.strip().split()[3]) # Displacement dim2
					Nodes.append(T_Displacement(dim1, dim2)) # Append displacements of the current node to the node list
					line = f.readline()
					
				if verbose: 
					print('Step {}: Dofs completed.\n'.format(nSteps))
					print('---------------------------------\n')
					
			# Read the stresses an strains at Gauss points
			elif line == 'Element output:\n': # String starts a list elements, GPs, strains and stresses
				line = f.readline()  # Cancel ---------- seperator
				line = f.readline()
				Elements = list() # Initialize/clear list of elements
				
				while line != '\n' and line !='\tR E A C T I O N S  O U T P U T:\n': # Strings that finish the list
#					nElement = line.strip().split()[2] # Element id
					line = f.readline()
					GPs = T_Element(list(), 0, 0, 0) # List of Gauss points

					while line != '\n' and line.strip().split()[0] == 'GP': # String that starts a new GP
#						nGP = int(line.strip().split()[1].split('.')[1]) # GP id
						tmp = [ float(i) for i in line.strip().split()[4:10] ] # Read the strains
						strain = T_Strains(tmp[0], tmp[1], tmp[2], tmp[3], tmp[4], tmp[5])
						line = f.readline()
						tmp = [ float(i) for i in line.strip().split()[1:7] ] # Read the stresses
						stress = T_Stresses(tmp[0], tmp[1], tmp[2], tmp[3], tmp[4], tmp[5])
						GPs.gp.append(T_GP(stress, strain)) # Append stresses and strains of the current GP to the GP list
						line = f.readline()
						
					Elements.append(GPs) # Append GP list of the current element to the element list 
				
				if verbose:
					print('Step {}: GPs completed.\n'.format(nSteps))
					print('---------------------------------\n')
			
				SimData.step.append(T_Step(Elements, Nodes)) # Append element and node list of the current step to the step list
			
			# Jump over the lines until we reach the next time step (Caught by if-clause)
			try:
				line = f.readline() # Will generate an error if files end is reached
			except:
				if verbose: print "End of file reached.\n"
				break # Break the 'while True' loop
	
	# -----------------------------------------------------------------------------
   
	# Averaging of strains and stress of GPs of each element
	for istep in range(len(SimData.step)):

		for ielement in range(len(SimData.step[istep].element)):
			# Initialization before each element
			stresses = np.array([0., 0., 0., 0., 0., 0.]) 
			strains = np.array([0., 0., 0., 0., 0., 0.])

			for igp in range(len(SimData.step[istep].element[ielement])):
				# Add up all data of all GPs
				stresses[:] += SimData.step[istep].element[ielement].gp[igp].stress[:]
				strains[:] += SimData.step[istep].element[ielement].gp[igp].strain[:]
			
			# Divide GP sum by number of GPs
			stresses /= len(SimData.step[istep].element[ielement])
			strains /= len(SimData.step[istep].element[ielement])
			# Replace the field (initialized with 0) with new information
			SimData.step[istep].element[ielement] = SimData.step[istep].element[ielement]._replace(avstress=T_Stresses(stresses[0],stresses[1],stresses[2],stresses[3],stresses[4],stresses[5]))
			SimData.step[istep].element[ielement] = SimData.step[istep].element[ielement]._replace(avstrain=T_Strains(strains[0],strains[1],strains[2],strains[3],strains[4],strains[5]))   
	return SimData

# -----------------------------------------------------------------------------

