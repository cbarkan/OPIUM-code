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

def valley_basis_find():
    r1 = np.array([1,-1*slope[0]/slope[1],0,0])
    coeff_vector = r1 + coeff_vector
    slope = slope_sign_find()*slope
    alpha = alpha_find()
    r1 = (coeff_vector - alpha*slope) - coeff_vector
    r1 = r1/np.linalg.norm(r1)
    
    r2 = np.array([1,0,-1*slope[0]/slope[2],0])
    coeff_vector = r2 + coeff_vector
    slope = slope_sign_find()*slope
    alpha = alpha_find()
    r2 = (coeff_vector - alpha*slope) - coeff_vector
    r2 = r1/np.linalg.norm(r2)
    
    r3 = np.array([1,0,0,-1*slope[0]/slope[3]])
    coeff_vector = r3 + coeff_vector
    slope = slope_sign_find()*slope
    alpha = alpha_find()
    r3 = (coeff_vector - alpha*slope) - coeff_vector
    r3 = r1/np.linalg.norm(r3)
    
    #Check for linear dependece
    dot12 = np.dot(r1,r2)
    dot13 = np.dot(r1,r3)
    dot23 = np.dot(r2,r3)
    if ( np.absolute(dot12) < 0.1 ) or ( np.absolute(dot13) < 0.1 ) or ( np.absolute(dot23) < 0.1 ):
        raise ValueError('Linearly dependent r vectors')
    
    R = np.array([r1,r2,r3])
    V = gs(R)
    A = V.T
    
    return A