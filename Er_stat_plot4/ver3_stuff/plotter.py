import numpy as np
import matplotlib.pylab as plt
from matplotlib.colors import LogNorm

mesh = np.load('zoom2_Er_stat.npy')
print mesh
"""
plt.figure(0)
plt.imshow(mesh, interpolation='none', origin = 'lower',extent=[-mesh_size,mesh_size,-mesh_size,mesh_size])
cb = plt.colorbar(orientation='horizontal')
cb.set_label('Er_stat')
plt.ylabel('Inner Box height (0.01 to 0.75 au)')
plt.xlabel('Outer Box height (0.75 to 1.5 au)')
"""
plt.figure(1)
plt.imshow(mesh, interpolation='none', origin = 'lower',extent=[0,4,-25,0], norm=LogNorm(vmin=1, vmax=5000))
cb = plt.colorbar(orientation='horizontal')
cb.set_label('Er_stat')
plt.ylabel('Inner Box height (0.01 to 0.75 au)')
plt.xlabel('Outer Box height (0.75 to 1.5 au)')
plt.show()
