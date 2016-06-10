
from matplotlib.colors import LogNorm
import numpy as np
import matplotlib.pyplot as plt


def make_single_inclusion(x,y,L):
    
    # allocate the matrix to hold the Green's function
    G = np.zeros([L,L]) 
    
    # iterate over all positions of the matrix
    for i in range(0,L): 
        for j in range(0,L):
            
            # calcualte the angle
            # if there is a problem with the arctan function, then use angle 0
            try: z = np.arctan(float(j-y)/float(i-x)) 
            except: z = 0. 
            
            # calculate distance between position and inclusion location
            r = np.sqrt((i-x)**2 + (j-y)**2) 
            
            # if distance is 0, do this to avoid dividing by 0
            if r>1e-12: G[i,j] += np.cos(4.*z)*np.power(r,-2.)
            else: G[i,j] += 0.
    
    return G


L = 128 

G = np.zeros([L,L]) 

# define x and y positions of all inclusion images
X = [-64,0,64,127,127+64]
Y = [-64,0,64,127,127+64]
    
# create array of images putting an inclusion at each position
for x in X:
    for y in Y:
        G += make_single_inclusion(x,y,L)
        
# put an inclusion at x_inc and y_inc and extract the stress from the array
# of images
x_inc = 20
y_inc = 50 
g = G[(L/2-x_inc):(L-x_inc), (L/2-y_inc):(L-y_inc)]

# create a figure
fig1 = plt.figure(figsize=(7,7),facecolor='white')
# add plot to the figure
ax1 = fig1.add_subplot(111)

# plot the matrix with log scale
im = ax1.imshow(np.abs(G), norm=LogNorm(vmin=1e-7, vmax=G.max()))

# add color bar
plt.colorbar(im)



#import greens_function
#from matplotlib.colors import LogNorm
#import numpy as np
#import matplotlib.pyplot as plt
#
## system size
#L = 64
#
## number of periodic images
#n = 1
#
## get the periodic array of Green's functions
#K = greens_function.get_periodic_kernel(L,n)
#
## get the stress field of a inclusion centered at (x,y)
#x = 2
#y = 25
#k = greens_function.get_inclusion_stress(K,x,y)
#
## the sum of the (internal) stress must be 0 (3rd Newton Law)
#print 'Sum of stress = ', k.sum()
#
## plot it
#fig1 = plt.figure(figsize=(7,7),facecolor='white')
#ax1 = fig1.add_subplot(111) 
#im = ax1.imshow(np.abs(k), norm=LogNorm(vmin=1e-7, vmax=K.max()))
#plt.colorbar(im)