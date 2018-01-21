#!/usr/bin/python
# coding: utf-8

import sys
import numpy as np
import mayavi
import math
from mayavi import mlab
# prepare some coordinates
voxels = np.genfromtxt(sys.argv[1], delimiter=',')
#dim = np.genfromtxt('voxels_dim.csv', delimiter=',')
#dim = dim.astype(np.int)
n = voxels.shape[0]
print n
dim = int(math.floor(np.power(n,1/3.0)))+1
print dim
voxels = np.reshape(voxels,(dim,dim,dim))
print voxels.shape
print (voxels>0).sum()

xx,yy,zz =  np.where(voxels>0)
s = voxels.flatten()
s = s[np.where(s>0)]
print xx.shape
print yy.shape
print zz.shape
print s.shape

# set the colors of each object
colors = np.empty(voxels.shape, dtype=object)
colors[voxels>1] = 'grey'
colors[voxels==1] = 'white'  
colors[voxels==8] = 'red'  
colors[voxels==7] = 'blue'  
colors[voxels==17] = 'green'  

voxels[voxels>0] = 1 
# and plot everything
mlab.points3d(xx, yy, zz,s,  
                     mode="cube",
	             colormap = "YlOrBr",
                     scale_mode='none',
                    scale_factor=1)

#mlab.contour3d(voxels,contours=2, transparent=True)
mlab.colorbar()
mlab.show()

