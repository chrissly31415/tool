import numpy as np
import mayavi
from mayavi import mlab
# prepare some coordinates
voxels = np.genfromtxt('voxels.csv', delimiter=',')
dim = np.genfromtxt('voxels_dim.csv', delimiter=',')
dim = dim.astype(np.int)
print voxels.shape
voxels = np.reshape(voxels,(dim[0],dim[1],dim[2]))
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

# and plot everything
#mlab.points3d(xx, yy, zz,s,  
#                     mode="cube",
#	             opacity=0.5,
#                     scale_factor=1)

voxels[voxels>0] = 1 
mlab.contour3d(voxels,contours=2, transparent=True)
mlab.show()

