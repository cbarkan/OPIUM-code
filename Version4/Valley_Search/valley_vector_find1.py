import numpy as np
import matplotlib.pylab as plt

def gs(X):
    #Returns matrix Y whose rows are an orthonormal basis of span{rows of X}
    Y = X.copy()
    Y[0] = Y[0]/np.linalg.norm(Y[0])
    for ii in range(1,Y.shape[0]):
        proj = 0
        for i in range(ii):
            proj += Y[i] * np.dot(Y[ii],Y[i])/np.dot(Y[i],Y[i])
        Y[ii] = Y[ii] - proj
        Y[ii] = Y[ii]/np.linalg.norm(Y[ii])
    return Y

def VtoC(x):
    return r[0] + np.dot(A,x)

def CtoV(x):
    #returns least square solution if x isn't in V-space
    v,residual = np.linalg.lstsq(A,x - r[0])[0:2]
    print 'IndtoV residual = %s' % residual
    return v

def CtoInd(coeff_vect):
    ind_vect = (coeff_vect - [lower_ext0,lower_ext1,lower_ext2,lower_ext3])/mesh_step
    return round(ind_vect[0]),round(ind_vect[1]),round(ind_vect[2]),round(ind_vect[3])

#_________________________________________________________________________

filename = 'zoom_n4-n8-n4-0'
mesh = np.load(filename + '_Eig1.npy') + np.load(filename + '_Norm1.npy')
c0_index = np.load(filename + '_c0_index.npy')
c1_index = np.load(filename + '_c1_index.npy')
c2_index = np.load(filename + '_c2_index.npy')
c3_index = np.load(filename + '_c3_index.npy')
lower_ext0 = c0_index[0]
upper_ext0 = c0_index[-1]
lower_ext1 = c1_index[0]
upper_ext1 = c1_index[-1]
lower_ext2 = c2_index[0]
upper_ext2 = c2_index[-1]
lower_ext3 = c3_index[0]
upper_ext3 = c3_index[-1]
mesh_step = 0.1

r_ind = np.zeros((6,4))
r = r_ind.copy()
r_ind[:,0:2] = np.array([[1,0],[2,1],[3,2],[4,4],[6,6],[7,8]])
for i in range(len(r_ind)):
    print np.min(mesh[r_ind[i,0],r_ind[i,1],:,:])
    a,b = np.where(mesh[r_ind[i,0],r_ind[i,1],:,:]==np.min(mesh[r_ind[i,0],r_ind[i,1],:,:]))
    r_ind[i,2:] = [a,b]
    r[i] = r_ind[i]*mesh_step + [lower_ext0,lower_ext1,lower_ext2,lower_ext3]

r_bar = r[1:4,:].copy()
for i in range(3):
    r_bar[i] = r_bar[i] - r[0]
    
v = gs(r_bar)
A = v.T




x = np.array([ 0.38729833,  0.31515762 , 0.22511258])
coeff_vect = VtoC(x)
print coeff_vect
a,b,c,d = CtoInd(coeff_vect)
print a,b,c,d
print mesh[a,b,c,d]

#print coeff_vect

"""
print CtoV(r[0])
print CtoV(r[1])
print CtoV(r[2])
print CtoV(r[3])
print CtoV(r[4])
print CtoV(r[5])
"""
"""
u = v
print np.dot(u[0],u[1])
print np.dot(u[2],u[1])
print np.dot(u[0],u[2])
print np.linalg.norm(u[0])
print np.linalg.norm(u[1])
print np.linalg.norm(u[2])
"""
#curve1 = mesh[3,2,:,5]
#surf = mesh[3,2,:,:]
#plt.imshow(surf,origin = 'lower')
#plt.plot(curve1)
#plt.show()