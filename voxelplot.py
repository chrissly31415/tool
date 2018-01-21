#!/usr/bin/python
# coding: utf-8

import matplotlib
import math
#matplotlib.use('agg')
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import sys
from matplotlib import cm, colors

print matplotlib.__version__

# prepare some coordinates
voxels = np.genfromtxt(sys.argv[1], delimiter=',')
n = voxels.shape[0]
dim = int(math.floor(np.power(n,1/3.0)))+1
print dim
voxels = np.reshape(voxels,(dim,dim,dim))

np.savetxt("tmp.csv",voxels,fmt='%5r')
print voxels.shape
print (voxels>0).sum()

# set the colors of each object
colors = np.empty(voxels.shape, dtype=object)
colors[voxels<1] = 'blue'
colors[voxels==1] = 'white'  
#colors[voxels==8] = 'red'  
#colors[voxels==7] = 'blue'  
#colors[voxels==17] = 'green'  
colors[voxels>1] = 'red'  

# and plot everything
fig = plt.figure()
ax = fig.gca(projection='3d')
#ax.voxels(voxels, facecolors=colors, edgecolor='k')

print cm.RdBu
ax.voxels(voxels, cmap=cm.RdBu, edgecolor='k')
maxlim = np.max(dim)
ax.set_xlim(0,maxlim)
ax.set_ylim(0,maxlim)
ax.set_zlim(0,maxlim)

plt.show()


