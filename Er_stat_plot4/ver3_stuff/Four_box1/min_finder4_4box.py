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

filename = 'n10-0-n5-5'
Total = np.load(filename + '_Total.npy')

box0_index = np.load(filename + '_box0_index.npy')
box1_index = np.load(filename + '_box1_index.npy')
box2_index = np.load(filename + '_box2_index.npy')
box3_index = np.load(filename + '_box3_index.npy')

print box0_index[0]

local_mins = np.zeros((1,5))
i = 0
size = Total.shape
for a in range(1,size[0]-1):
    for b in range(1,size[1]-1):
        for c in range(1,size[2]-1):
            for d in range(1,size[3]-1):
                err = Total[a,b,c,d]
                err_neighbor = min_directNeighbor_find(a,b,c,d)
                if err < err_neighbor:
                    i0 = box0_index[a]
                    i1 = box1_index[b]
                    i2 = box2_index[c]
                    i3 = box3_index[d]
                    local_mins[-1,:] = [i0,i1,i2,i3,err]
                    local_mins = np.append(local_mins,np.zeros((1,5)),axis=0)
                    
np.set_printoptions(suppress=True)
print local_mins
