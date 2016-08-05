import numpy as np
import matplotlib.pylab as plt
from matplotlib.colors import LogNorm

filename = 'n4_n10-n6-n6-n2'

mesh1 = np.load(str(filename) + '_Eig.npy')

mesh2 = np.load(str(filename) + '_Norm.npy')




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
min_value = 1
max_value = 1000
#Log
plt.figure(1)
plt.imshow(mesh1, interpolation='none', origin = 'lower',extent=[lower_ext1,upper_ext1,lower_ext0,upper_ext0], norm=LogNorm(vmin=min_value, vmax=max_value))
cb = plt.colorbar(orientation='horizontal')
cb.set_label('Er_stat')
plt.ylabel('Inner Box height (0.01 to 0.75 au)')
plt.xlabel('Outer Box height (0.75 to 1.5 au)')

plt.figure(2)
plt.imshow(mesh2, interpolation='none', origin = 'lower',extent=[lower_ext1,upper_ext1,lower_ext0,upper_ext0], norm=LogNorm(vmin=min_value, vmax=max_value))
cb = plt.colorbar(orientation='horizontal')
cb.set_label('Er_stat')
plt.ylabel('Inner Box height (0.01 to 0.75 au)')
plt.xlabel('Outer Box height (0.75 to 1.5 au)')

plt.show()
