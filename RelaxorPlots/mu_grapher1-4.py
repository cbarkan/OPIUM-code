import numpy as np

filenames = ['mu_total1to250.npy','mu_total251to475.npy','mu_total476to629.npy']

mu = np.load(filenames[0])
for i in range(1,len(filenames)):
    mu = np.vstack((mu,np.load(filenames[i])))

"""                              
mu1 = np.load(filename1)
mu2 = np.load(filename2)
mu3 = np.load(filename2)
mu = np.vstack((mu1,mu2))
"""

x = mu[:,0]
y = mu[:,1]
z = mu[:,2]
print np.min(x)
print np.max(x)
print np.min(y)
print np.max(y)
print np.min(z)
print np.max(z)


import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.plot(x,y,z)
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('z')
ax.text(x[0],y[0],z[0],'start')
ax.scatter(x[0],y[0],z[0],s=40,c='red')
ax.text(x[-1],y[-1],z[-1],'end')
ax.scatter(x[-1],y[-1],z[-1],s=40,c='red')

#Find dump:
dump_num_list = [245]

for dump_num in dump_num_list:
    label = 'dump %s' % dump_num
    ax.text(x[dump_num],y[dump_num],z[dump_num],label)
    ax.scatter(x[dump_num],y[dump_num],z[dump_num],s=40,c='red')

plt.show()