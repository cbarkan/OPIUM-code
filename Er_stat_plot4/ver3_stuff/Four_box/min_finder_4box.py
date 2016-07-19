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


Eig = np.load('small_test1_Eig.npy')
Norm = np.load('small_test1_Norm.npy')
Total = Eig + Norm


local_mins = np.zeros((1,5))
i = 0
size = Total.shape
for a in range(1,size[0]-1):
    for b in range(1,size[1]-1):
        for c in range(1,size[2]-1):
            for d in range(1,size[3]-1):
                err = Total[a,b,c,d]
                err_neighbor = min_directNeighbor_find(a,b,c,d)
                print'%s %s %s %s   err = %s' % (a, b, c, d, err)
                if err < err_neighbor:
                    local_mins[-1,:] = [a,b,c,d,err]
                    local_mins = np.append(local_mins,np.zeros((1,5)),axis=0)
                    
np.set_printoptions(suppress=True)
print local_mins