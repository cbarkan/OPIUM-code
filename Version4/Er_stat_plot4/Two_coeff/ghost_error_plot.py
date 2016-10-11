import numpy as np
import matplotlib.pylab as plt

filename = 'test1'



box0_index = np.load(filename + '_c0_index.npy')
box1_index = np.load(filename + '_c1_index.npy')
lower_ext0 = box0_index[0]
upper_ext0 = box0_index[-1]
lower_ext1 = box1_index[0]
upper_ext1 = box1_index[-1]
Err = np.load(filename+'_Eig1.npy')
l1,l2 = Err.shape

mesh = Err.copy()

for ii in range(l1):
	for i in range(l2):
		if Err[ii,i] == np.inf: #error
			mesh[ii,i] = 0
		elif Err[ii,i] == -1*np.inf: #ghost
			mesh[ii,i] = 1
		else:
			mesh[ii,i] = np.inf

plt.imshow(mesh, interpolation='none', origin = 'lower',extent=[lower_ext1,upper_ext1,lower_ext0,upper_ext0],vmin=-0.5, vmax=1.5)
plt.ylabel('a0')
plt.xlabel('a1')
plt.title('V')
plt.show()