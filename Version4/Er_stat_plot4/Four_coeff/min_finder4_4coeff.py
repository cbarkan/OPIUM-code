import numpy as np

def min_directNeighbor_find(a,b,c,d):
    #This does not look for box settings diagonal to [a,b,c,d]
    #kb_find3 doesn't search for diagonals either.
    #I should search for minima that are diagonal to the minima found
    a_plus = Total[a+1,b,c,d]
    a_minus = Total[a-1,b,c,d]
    b_plus = Total[a,b+1,c,d]
    b_minus = Total[a,b-1,c,d]
    c_plus = Total[a,b,c+1,d]
    c_minus = Total[a,b,c-1,d]
    d_plus = Total[a,b,c,d+1]
    d_minus = Total[a,b,c,d-1]
    return min(a_plus,a_minus,b_plus,b_minus,c_plus,c_minus,d_plus,d_minus)

filename = 'zoom_n4-n8-n4-0'
Total = np.load(filename + '_Eig1.npy') + np.load(filename + '_Norm1.npy')

c0_index = np.load(filename + '_c0_index.npy')
c1_index = np.load(filename + '_c1_index.npy')
c2_index = np.load(filename + '_c2_index.npy')
c3_index = np.load(filename + '_c3_index.npy')

local_mins = np.zeros((1,5))
i = 0
size = Total.shape
for a in range(1,size[0]-1):
    for b in range(1,size[1]-1):
        for c in range(1,size[2]-1):
            for d in range(1,size[3]-1):
                err = Total[a,b,c,d]
                err_neighbor = min_directNeighbor_find(a,b,c,d)
                if err <= err_neighbor:
                    i0 = c0_index[a]
                    i1 = c1_index[b]
                    i2 = c2_index[c]
                    i3 = c3_index[d]
                    local_mins[-1,:] = [i0,i1,i2,i3,err]
                    local_mins = np.append(local_mins,np.zeros((1,5)),axis=0)
                    
np.set_printoptions(suppress=True)
print local_mins




abs_min = np.min(Total)
print 'absolute min = %s' % abs_min
a,b,c,d = np.where(Total[:,:,:,:] == abs_min)[0:4]
c0_min = c0_index[a[0]]
c1_min = c1_index[b[0]]
c2_min = c2_index[c[0]]
c3_min = c3_index[d[0]]
print c0_min, c1_min, c2_min, c3_min