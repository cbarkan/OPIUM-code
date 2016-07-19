import numpy as np

a = np.load('Eig_mesh.npy')
b = np.load('Norm_mesh.npy')
err = a + b
err5by5 = err #[15:25,15:25]

minimum = np.min(err5by5)
(l,r) = np.where(err5by5 == minimum)
print l
print r
print minimum