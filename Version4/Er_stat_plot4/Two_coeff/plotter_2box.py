import numpy as np
import matplotlib.pylab as plt
from matplotlib.colors import LogNorm

filename = 'cu'

mesh1 = np.load(str(filename) + '_Eig5.npy') + np.load(str(filename) + '_Norm5.npy')

"""
mesh2 = np.load(str(filename) + '_Eig2.npy')
mesh3 = np.load(str(filename) + '_Eig3.npy')
mesh4 = np.load(str(filename) + '_Eig4.npy')
mesh5 = np.load(str(filename) + '_Eig5.npy')
"""

#mesh1 = np.load(str(filename) + '_Norm1.npy')
"""
mesh2 = np.load(str(filename) + '_Norm2.npy')
mesh3 = np.load(str(filename) + '_Norm3.npy')
mesh4 = np.load(str(filename) + '_Norm4.npy')
mesh5 = np.load(str(filename) + '_Norm5.npy')
"""


box0_index = np.load(str(filename) + '_c0_index.npy')
box1_index = np.load(str(filename) + '_c1_index.npy')
lower_ext0 = box0_index[0]
upper_ext0 = box0_index[-1]
lower_ext1 = box1_index[0]
upper_ext1 = box1_index[-1]

"""
###Zoom: Comment block out if not zooming
lower_ext0 = 0.24
upper_ext0 = 0.25
lower_ext1 = 4.055
upper_ext1 = 4.065

a = np.where(box0_index == lower_ext0)
print a

lower_ext0_index = np.int(np.where(abs(box0_index-lower_ext0)<0.00000001)[0])
upper_ext0_index = np.int(np.where(abs(box0_index-upper_ext0)<0.00000001)[0])
lower_ext1_index = np.int(np.where(abs(box1_index-lower_ext1)<0.00000001)[0])
upper_ext1_index = np.int(np.where(abs(box1_index-upper_ext1)<0.00000001)[0])

mesh = mesh[lower_ext0_index:upper_ext0_index,lower_ext1_index:upper_ext1_index]
###End zoom block
"""
"""
#Lin
plt.figure(0)
plt.imshow(mesh, interpolation='none', origin = 'lower',extent=[lower_ext1,upper_ext1,lower_ext0,upper_ext0])
cb = plt.colorbar(orientation='horizontal')
cb.set_label('Er_stat')
plt.ylabel('Inner Box height (0.01 to 0.75 au)')
plt.xlabel('Outer Box height (0.75 to 1.5 au)')
"""
title = filename[0:2]
if title[1] == '_':
    title = title[0]

title = 'Cu'
xlab = 'a0'
ylab = 'a1'

#Log
plt.figure(1)
plt.imshow(mesh1, interpolation='none', origin = 'lower',extent=[lower_ext1,upper_ext1,lower_ext0,upper_ext0], norm=LogNorm(vmin=1, vmax=150))
cb = plt.colorbar(orientation='horizontal')
cb.set_label('Er_stat')
plt.ylabel(xlab)
plt.xlabel(ylab)
plt.title(title)

"""
plt.figure(2)
plt.imshow(mesh2, interpolation='none', origin = 'lower',extent=[lower_ext1,upper_ext1,lower_ext0,upper_ext0], norm=LogNorm(vmin=1, vmax=5000))
cb = plt.colorbar(orientation='horizontal')
cb.set_label('Er_stat')
plt.ylabel(xlab)
plt.xlabel(ylab)
plt.title(title)

plt.figure(3)
plt.imshow(mesh3, interpolation='none', origin = 'lower',extent=[lower_ext1,upper_ext1,lower_ext0,upper_ext0], norm=LogNorm(vmin=1, vmax=5000))
cb = plt.colorbar(orientation='horizontal')
cb.set_label('Er_stat')
plt.ylabel(xlab)
plt.xlabel(ylab)
plt.title(title)

plt.figure(4)
plt.imshow(mesh4, interpolation='none', origin = 'lower',extent=[lower_ext1,upper_ext1,lower_ext0,upper_ext0], norm=LogNorm(vmin=1, vmax=5000))
cb = plt.colorbar(orientation='horizontal')
cb.set_label('Er_stat')
plt.ylabel(xlab)
plt.xlabel(ylab)
plt.title(title)

plt.figure(5)
plt.imshow(mesh5, interpolation='none', origin = 'lower',extent=[lower_ext1,upper_ext1,lower_ext0,upper_ext0], norm=LogNorm(vmin=1, vmax=5000))
cb = plt.colorbar(orientation='horizontal')
cb.set_label('Er_stat')
plt.ylabel(xlab)
plt.xlabel(ylab)
plt.title(title)


#Ghosts
mesh_ghost = np.isneginf(mesh1)
plt.figure(6)
plt.imshow(mesh_ghost, interpolation='none', origin = 'lower',extent=[lower_ext1,upper_ext1,lower_ext0,upper_ext0], norm=LogNorm(vmin=1, vmax=5000))
cb = plt.colorbar(orientation='horizontal')
cb.set_label('Er_stat')
plt.ylabel(xlab)
plt.xlabel(ylab)
plt.title(title + 'ghost')

#OPIUM error
mesh_err = np.isneginf(-1*mesh1)
plt.figure(7)
plt.imshow(mesh_err, interpolation='none', origin = 'lower',extent=[lower_ext1,upper_ext1,lower_ext0,upper_ext0], norm=LogNorm(vmin=1, vmax=5000))
cb = plt.colorbar(orientation='horizontal')
cb.set_label('Er_stat')
plt.ylabel(xlab)
plt.xlabel(ylab)
plt.title(title + 'err')
"""

plt.show()
