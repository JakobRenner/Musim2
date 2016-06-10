import numpy as np


def get_periodic_kernel(L,nIm):

    stressArray = np.zeros([2*L+1,2*L+1])
    dx = np.zeros(4)
    dy = np.zeros(4)
    b = np.zeros(4)
    
    for n in range(2*nIm+1):
        for i in range(2*L+1):
            for j in range(2*L+1):
                
                dx[0] = i-(L+0.5)+L*(nIm+1-n)
                dy[0] = j-L
                b[0] = 1

                dy[1] = i-L                      
                dx[1] = j-(L+0.5)+L*(nIm+1-n)    
                b[1] = 1

                dx[2] = i-(L-0.5)+L*(nIm+1-n)
                dy[2] = j-L
                b[2] = -1

                dy[3] = i-L                      
                dx[3] = j-(L-0.5)+L*(nIm+1-n)    
                b[3] = -1
                
                denom = np.power(np.cosh(2.*3.14159*dx/L)-np.cos(2.*3.14159*dy/L),2.)
                addStress = b*2.*3.14159*dx/L*(np.cosh(2.*3.14159*dx/L)*np.cos(2.*3.14159*dy/L)-1)/denom

                stressArray[i,j] = stressArray[i,j] + sum(addStress)
                
    norm = -1/stressArray.min()
    stressArray = stressArray*norm

    return stressArray

def get_inclusion_stress(K,x,y):
    L = int(K.shape[0]/2)
    return K[(L-x):(2*L-x), (L-y):(2*L-y)]
