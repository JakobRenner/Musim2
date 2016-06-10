import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

#theta = np.arange(0, 2*np.pi, res)

n = 64
x = 32
y = 32

G = np.zeros([n, n])

for i in range( n ):
    for j in range( n ):
        try:
            theta = np.arctan( (j-y) / (i-x) )
        except: theta = 0.
        r = np.sqrt( (i - x)**2 + (j - y)**2 )


        if r < 1e-12:
            G[i][j] = 0
        else:
            G[i][j] = np.cos(4 * theta) / r**2

fig1 = plt.figure(figsize =(7,7), facecolor = 'white')
ax1 = fig1.add_subplot(111)
im = ax1.matshow( np.abs(G), norm=LogNorm(vmin=1e-7, vmax = G.max()), cmap = 'winter' )
plt.colorbar(im)
plt.show()
