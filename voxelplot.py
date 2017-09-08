import matplotlib
#matplotlib.use('agg')
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D

print matplotlib.__version__

# prepare some coordinates
voxels = np.genfromtxt('voxels.csv', delimiter=',')
dim = np.genfromtxt('voxels_dim.csv', delimiter=',')
dim = dim.astype(np.int)
print dim
voxels = np.reshape(voxels,(dim[0],dim[1],dim[2]))
np.savetxt("tmp.csv",voxels,fmt='%5r')
print voxels.shape
print (voxels>0).sum()

# set the colors of each object
colors = np.empty(voxels.shape, dtype=object)
colors[voxels>1] = 'grey'
colors[voxels==1] = 'white'  
colors[voxels==8] = 'red'  
colors[voxels==7] = 'blue'  
colors[voxels==17] = 'green'  

# and plot everything
fig = plt.figure()
ax = fig.gca(projection='3d')
#ax.voxels(voxels, facecolors=colors, edgecolor='k')
ax.voxels(voxels, facecolors=colors, edgecolor='k')
maxlim = np.max(dim)
ax.set_xlim(0,maxlim)
ax.set_ylim(0,maxlim)
ax.set_zlim(0,maxlim)

plt.show()


