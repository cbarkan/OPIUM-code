import numpy as np
import matplotlib.pylab as plt

filename = 'test1'
mesh = np.load(filename+'_Eig.npy') + np.load(filename+'_Norm.npy')
c0_index = np.load(filename+'_c0_index.npy')
c1_index = np.load(filename+'_c1_index.npy')

section = mesh[:,22]
print section
plt.plot(section)
plt.show()