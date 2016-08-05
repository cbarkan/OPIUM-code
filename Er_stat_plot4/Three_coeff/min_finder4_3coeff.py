import numpy as np

def min_directNeighbor_find(a,b,c):
    #This does not look for coeffs diagonal to [a,b,c]
    a_plus = Total[a+1,b,c]
    a_minus = Total[a-1,b,c]
    b_plus = Total[a,b+1,c]
    b_minus = Total[a,b-1,c]
    c_plus = Total[a,b,c+1]
    c_minus = Total[a,b,c-1]
    return min(a_plus,a_minus,b_plus,b_minus,c_plus,c_minus)

filename = 'test1'
Total = np.load(filename + '_Eig1.npy') + np.load(filename + '_Norm1.npy')

c0_index = np.load(filename + '_c0_index.npy')
c1_index = np.load(filename + '_c1_index.npy')
c2_index = np.load(filename + '_c2_index.npy')

#print cx0_index[0]

local_mins = np.zeros((1,4))
i = 0
size = Total.shape
for a in range(1,size[0]-1):
    for b in range(1,size[1]-1):
        for c in range(1,size[2]-1):
            err = Total[a,b,c]
            err_neighbor = min_directNeighbor_find(a,b,c)
            if err < err_neighbor:
                i0 = c0_index[a]
                i1 = c1_index[b]
                i2 = c2_index[c]
                local_mins[-1,:] = [i0,i1,i2,err]
                local_mins = np.append(local_mins,np.zeros((1,4)),axis=0)
                    
np.set_printoptions(suppress=True)
np.savetxt(filename + '_Eig1_mins.txt',local_mins,'%5.3f',delimiter='    ')
print local_mins
print np.min(local_mins[0:-1,3])