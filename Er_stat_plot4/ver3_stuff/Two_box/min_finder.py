import numpy as np

mesh = np.load('Er_stat.npy')
Left_index = np.load('Left_index.npy')
Right_index = np.load('Right_index.npy')
(l,r) = np.where(mesh == np.min(mesh))
l_height = Left_index[l]
r_height = Right_index[r]
print l_height
print r_height
