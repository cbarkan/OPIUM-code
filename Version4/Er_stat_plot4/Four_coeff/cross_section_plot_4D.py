import numpy as np
import matplotlib.pylab as plt

filename = 'superzoom1_n3p8-n7p6-n3p8-np2'
mesh = np.load(filename+'_Eig1.npy') + np.load(filename+'_Norm1.npy')
#c0_index = np.load(filename+'_c0_index.npy')
#c1_index = np.load(filename+'_c1_index.npy')

section = mesh[5,5,:,3]
print section
plt.plot(section,linestyle='None',marker='o')
plt.show()