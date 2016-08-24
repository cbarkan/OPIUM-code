import numpy as np
import matplotlib.pylab as plt

filename = 'trial1'
mesh = np.load(filename+'_Eig1.npy') + np.load(filename+'_Norm1.npy')
#c0_index = np.load(filename+'_c0_index.npy')
#c1_index = np.load(filename+'_c1_index.npy')

section = mesh[10,10,:]
print section
plt.plot(section)
plt.show()