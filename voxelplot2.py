#!/usr/bin/python
# coding: utf-8

import sys
import numpy as np
import mayavi
import math
from mayavi import mlab

import matplotlib.pyplot as plt

def analyze_voxels(voxels):
   print "voxel.shape:",voxels.shape
   print "voxel sum>0.0:",(voxels>0).sum()
   print "voxel sum<0.0:",(voxels<0).sum()
   print "voxel sum=0.0:",(voxels==0).sum()
   print "voxel max:",voxels.max()
   print "voxel min:",voxels.min()

# prepare some coordinates
voxels = np.genfromtxt(sys.argv[1], delimiter=',')
#dim = np.genfromtxt('voxels_dim.csv', delimiter=',')
#dim = dim.astype(np.int)
n = voxels.shape[0]
print "shape:",n
dim = int(math.floor(np.power(n,1/3.0)))+1
print "dimension:",dim

voxels = np.reshape(voxels,(dim,dim,dim))

#extract isosurfaces

if len(sys.argv)>2:
 	vmax= float(sys.argv[3]) 
 	vmin= float(sys.argv[2])
else:
	vmax = 20.0  
	vmin = 0.5

print("vmax",vmax)
print("vmin",vmin)

analyze_voxels(voxels) 

plt.hist(voxels.flatten(),bins=30, alpha=0.5)

#vmax = voxels.flatten().max()
idx1 = np.logical_and(voxels<vmax,voxels>vmin)
idx2 = np.logical_or(voxels<vmin,voxels>vmax)
voxels[idx2] = 0.0
plt.hist(voxels.flatten(),bins=30, alpha=0.5)
plt.show()

analyze_voxels(voxels) 

lb=0.0
xx,yy,zz =  np.where(voxels!=lb)
s = voxels.flatten()
s = s[np.where(s!=lb)]
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

# and plot everything
mlab.points3d(xx, yy, zz,s,  
                     mode="cube",
#                     vmax=200,
#	             vmin=0.0,
	             colormap = "YlOrBr",
                     scale_mode='none',
                    scale_factor=1)

#mlab.contour3d(voxels,contours=2, transparent=True)
mlab.colorbar()
mlab.show()

