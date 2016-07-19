import numpy as np
mesh = np.load('test1_Total.npy')
c0_index = np.load('test1_c0_index.npy')
c1_index = np.load('test1_c1_index.npy')


coeff0 = 11
coeff1 = 11

c0 = np.where(c0_index == coeff0)
c1 = np.where(c1_index == coeff1)

Er = mesh[c0,c1]
print Er