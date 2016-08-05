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
    print a_plus
    print a_minus
    print b_plus
    print b_minus
    print c_plus
    print c_minus
    print d_plus
    print d_minus
    """
    if a_plus == min(a_plus,a_minus,b_plus,b_minus,c_plus,c_minus,d_plus,d_minus):
        return np.array([a+1,b,c,d]), a_plus
        """
    return min(a_plus,a_minus,b_plus,b_minus,c_plus,c_minus,d_plus,d_minus)


filename = 'zoom_n4-n8-n4-0'
Total = np.load(filename + '_Eig1.npy') + np.load(filename + '_Norm1.npy')

c0_index = np.load(filename + '_c0_index.npy')
c1_index = np.load(filename + '_c1_index.npy')
c2_index = np.load(filename + '_c2_index.npy')
c3_index = np.load(filename + '_c3_index.npy') 

c0,c1,c2,c3 = -4.0,-8.0,-4.0,0.0
a = np.where(c0_index == c0)[0]
b = np.where(c1_index == c1)[0]
c = np.where(c2_index == c2)[0]
d = np.where(c3_index == c3)[0]

print 'Current Error = %s' % Total[a,b,c,d]
new_min = min_directNeighbor_find(a,b,c,d)
print new_min