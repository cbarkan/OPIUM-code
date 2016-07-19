import numpy as np

mesh = np.load('test1_Total.npy')
c0_index = np.load('test1_c0_index.npy')
c1_index = np.load('test1_c1_index.npy')
minimum = np.min(mesh)
print minimum
(c0,c1) = np.where(mesh == minimum)
c0_height = c0_index[c0]
c1_height = c1_index[c1]
print c0_height
print c1_height
