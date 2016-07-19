import numpy as np
import matplotlib.pylab as plt
from matplotlib.colors import LogNorm

mesh = np.load('n.15-n.05-n.02-n.12_Total.npy')
box0_index = np.load('n.15-n.05-n.02-n.12_Left_index.npy')
box1_index = np.load('n.15-n.05-n.02-n.12_Right_index.npy')
lower_ext0 = box0_index[0]
upper_ext0 = box0_index[-1]
lower_ext1 = box1_index[0]
upper_ext1 = box1_index[-1]

plt.figure(0)
plt.imshow(mesh, interpolation='none', origin = 'lower',extent=[lower_ext0,upper_ext0,lower_ext1,upper_ext1])
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
"""

plt.show()
