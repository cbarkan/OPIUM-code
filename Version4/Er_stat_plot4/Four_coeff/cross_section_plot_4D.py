import numpy as np
import matplotlib.pylab as plt

filename = 'step3'
mesh = np.load(filename+'_Eig1.npy') + np.load(filename+'_Norm1.npy')
#c0_index = np.load(filename+'_c0_index.npy')
#c1_index = np.load(filename+'_c1_index.npy')

#1D Cross section
"""
section = mesh[0,0,0:30,7]
print section
plt.plot(section,linestyle='None',marker='o')
c0_index = np.load(filename + '_c0_index.npy')
c1_index = np.load(filename + '_c1_index.npy')
c2_index = np.load(filename + '_c2_index.npy')
c3_index = np.load(filename + '_c3_index.npy')
title = '%s, %s, \n %s to %s, \n %s' % (c0_index[0],c1_index[0],c2_index[0],c2_index[-1],c3_index[5])
plt.title(title)
plt.show()
"""

#2D cross section

c2_lowerBound = 15
c2_upperBound = 25
c3_lowerBound = 15
c3_upperBound = 25
section = mesh[0,0,c2_lowerBound:c2_upperBound,c3_lowerBound:c3_upperBound]

c0_index = np.load(filename + '_c0_index.npy')
c1_index = np.load(filename + '_c1_index.npy')
c2_index = np.load(filename + '_c2_index.npy')
c3_index = np.load(filename + '_c3_index.npy')

lower_ext0 = c0_index[0]
upper_ext0 = c0_index[-1]
lower_ext1 = c1_index[0]
upper_ext1 = c1_index[-1]
lower_ext2 = c2_index[c2_lowerBound]
upper_ext2 = c2_index[c2_upperBound]
lower_ext3 = c3_index[c3_lowerBound]
upper_ext3 = c3_index[c3_upperBound]

l = plt.imshow(section,interpolation='none', origin = 'lower',extent=[lower_ext3,upper_ext3,lower_ext2,upper_ext2])
#cb = plt.colorbar(l,orientation='horizontal')
plt.ylabel('Coefficient 2')
plt.xlabel('Coefficient 3')
plt.show()