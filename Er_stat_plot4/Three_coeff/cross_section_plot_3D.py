import numpy as np
import matplotlib.pylab as plt

filename = 'test1'
mesh = np.load(filename+'_Eig1.npy') + np.load(filename+'_Norm1.npy')
c0_index = np.load(filename+'_c0_index.npy')
c1_index = np.load(filename+'_c1_index.npy')

section = mesh[9,6,:]
print section
plt.plot(section)
plt.show()