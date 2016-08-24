
import numpy as np
import matplotlib.pylab as plt
from subprocess import Popen
import threading


###Inputs###
filename = 'v' #(must be string) Write this as string without a file extension (i.e. 'Fe' is correct, 'Fe.param' is incorrect)
z = 23 #Atomic number of element
tconfig_num = 5 #(must be integer) Number of test configurations, this number must match the corresponding number in the .param file
alpha_step = 5.0 #(must be float, if integer is desired, write 2.0, for example) Initial step size in the search for optimal alpha
local_orb = 's' #(must be string, either 's', 'p', 'd', or 'f') local orbital setting for augmentation function. 
a = 0.0001 #Parameter which defines grid point spacing. Do not change unless you make a corresponding change in OPIUM
b = 0.013 #Parameter which defines grid point spacing. Do not change unless you make a corresponding change in OPIUM
timeout_sec = 30

alpha_step = 100000.0
slopeCheckAlpha = 0.0000005
resolution = 50.0 #boxes per a.u.
x_min = 0.01
x_max = 1.5
############

#Initialize grid
grid_size = int((x_max - x_min)*resolution)
grid = np.zeros((grid_size+1,2))
grid[0,0] = np.rint(np.log(x_min/(a*(z**(-1.0/3.0))))/b + 1)
for i in np.arange(1,grid_size+1): #Set box bounds
    grid[i,0] = i*(x_max-x_min)/grid_size #UNITS: a.u.
    grid[i,0] = np.log(grid[i,0]/(a*(z**(-1.0/3.0))))/b + 1 #Converts to g.p units
    grid[i,0] = np.rint(grid[i,0]) #Ensures integer values of g.p


def fourier_cos(x,coeff_vector): #Evaluates F(x) where F is a fourier cosine series with the given coefficients
    f = 0.0
    L = x_max
    for i in range(coeff_terms):
        f += coeff_vector[i]*np.cos(x*i*np.pi/L)
    return f

def fit_to_func(coeff_vector_test,form):
    if form=='fourier_cos':
        for i in range(grid_size):
            grid[i,1] = fourier_cos(x_min + (i+0.5)/resolution,coeff_vector_test)
    elif form=='poly':
        for i in range(grid_size):
            grid[i,1] = poly(x_min + (i+0.5)/resolution,coeff_vector_test)
    else:
        return 'Error: functional form not properly specified'
    return grid

def setKB(grid): #Sets [KBdesign] in .param file to specified grid
    #Get param file text
    param_file = open(filename + '.param', 'r')
    param = param_file.readlines()
    param_file.close

    #Open param file for rewrite
    param_file = open(filename + '.param', 'w')

    index = 0
    for jj in range(len(param) - 2):
        line = param[index]
        if line.startswith('[K'):
            param_file.write('[KBdesign] \n' + local_orb + '\n' + str(grid_size) + '\n')
            subindex = 0
            for i in np.arange(grid_size): #This adds all the new boxes
                param_file.write("gp " + str(grid[subindex,0]) + " " + str(grid[subindex+1,0] - 1) + " " + str(grid[subindex,1]) + "\n")
                subindex += 1
            
            index +=3 #jumps down to pre-existing boxes
            
        elif line.startswith("gp") or line.startswith("au"): #This deletes all previus boxes
            index += 1
        else: #This copies the rest of the file as is
            param_file.write(line)
            index += 1
                
    param_file.close
    return

def VtoC(x):
    return origin_displacement + np.dot(A,x)

def err_check(alpha,coord_vector,slope,space_type): #Runs OPIUM, checks test config errors and returns their sum

    if space_type == 'coeff':
        test_vector = coord_vector.copy()
        test_vector = test_vector - alpha*slope
        grid = fit_to_func(test_vector,'fourier_cos')
        setKB(grid)
    elif space_type == 'valley':
        test_vector = VtoC(coord_vector - alpha*slope)
        grid = fit_to_func(test_vector,'fourier_cos')
        setKB(grid)
    else:
        raise ValueError('space_type must be "coeff" or "valley"')
    
    #Use threading to cancel OPIUM if it takes too long
    proc = Popen('./opium ' + filename + ' ' + filename + '.log all rpt', shell=True)
    proc_thread = threading.Thread(target=proc.communicate)
    proc_thread.start()
    proc_thread.join(timeout_sec)
    if proc_thread.is_alive():
        # Process still running - kill it and raise timeout error
        try:
            proc.kill()
        except OSError, e:
            # The process finished between the `is_alive()` and `kill()`
            return proc.returncode

    rpt_file = open(filename + '.rpt', 'a+')
    rpt = rpt_file.readlines()
    if rpt[-1] == '#kb_find read-receipt': #Error occured in OPIUM and no new .rpt file was produced
        Er_stat = float('inf') #Returns inf so that current KB box is not chosen
    else: #No error occured
        Er = np.ones((tconfig_num,2))
        Er_index = 0
        DD = np.ones((tconfig_num*(tconfig_num+1)/2,1))
        DD_index = 0
        for line in rpt:
            if line.startswith("  !!ERROR!!"):
                #Ghosts in psp
                return float('inf')
            if line.startswith(" AE-NL-  t"):
                Er[Er_index,0] = float(line[22:38])
                Er[Er_index,1] = float(line[42:58])
                Er_index += 1
            elif line.startswith(" AE-NL-   ") and (not line.startswith(" AE-NL-   i")):
                DD[DD_index,0] = abs(float(line[22:33]))
                DD_index += 1
        
        rpt_file.write('#kb_find read-receipt')
        rpt_file.close

        #Er_stat = sum(Er[:,0]) #Sum of test config eigenvalue errors
        #Er_stat = sum(sum(DD)) #Sum of energy differences
        #Er_stat = sum(sum(Er)) #Sum of eigenvalue errors and norm errors
        #Er_stat = sum(sum(Er))/tconfig_num + sum(sum(DD))/(tconfig_num*(tconfig_num-1)/2) #Sum of avg eigenvalue error, avg norm error, and avg energy difference
        Er_stat = Er[0,0] + Er[0,1]
        
        global best_coeff_vector
        global best_coeff_vector_error
        global best_coeff_iter
        if Er_stat < best_coeff_vector_error:
            best_coeff_vector_error = Er_stat
            best_coeff_vector = test_vector
            best_coeff_iter = prim_iter
            print 'New best coeff_vector!!'
            
    return Er_stat


phi = 1.6180339887
def alpha_find(coord_vector,alpha_step,space_type):
    #aL / ErL == alpha_lower / Error_lower
    #amL / ErmL == alpha_middle-lower / Error_middle-lower
    #amU / ErmU == alpha_middle-upper / Error_middle-upper
    #aL / ErU == alpha_upper / Error_upper
    """
    if space_type == 'coeff':
        coord_vector = coeff_vector
    elif space_type == 'valley':
        coord_vector = valley_vector
    else:
        raise ValueError('space_type must be "coeff" or "valley"')
    """

    slope_length = np.linalg.norm(slope)
    if slope_length==0:
        print 'Exiting alpha_find: slope_length == 0. Returning alpha = 0'
        return 0
    aL = 0
    amL = 0
    ErL = err_check(aL,coord_vector,slope,space_type)
    ErmL = ErL
    print 'testing1:'
    print coord_vector
    print ErL
    
    #Let alpha travel
    aU = alpha_step
    ErU = err_check(aU,coord_vector,slope,space_type)
    print ErU
    
    if ErU < ErmL:
        while ErU < ErmL:
            ErL = ErmL
            aL = amL
            ErmL = ErU
            amL = aU
            aU += alpha_step
            ErU = err_check(aU,coord_vector,slope,space_type)
            print ErU
    else:
        print 'Alpha didnt travel, reducing alpha_step'
        for i in range(10000):
            aU = aU/4
            alpha_step = alpha_step/4
            print 'alpha_step = %s' % alpha_step
            ErU = err_check(aU,coord_vector,slope,space_type)
            print ErU
            if ErU < ErmL:
                print ErU
                while ErU < ErmL:
                    ErL = ErmL
                    aL = amL
                    ErmL = ErU
                    amL = aU
                    aU += alpha_step
                    ErU = err_check(aU,coord_vector,slope,space_type)
                    print ErU
                break
            elif alpha_step < 0.000000001*slopeCheckAlpha/slope_length:
                print 'No alpha found: returning alpha = %s' % 0
                return 0 #If loop completes, no good alpha found, return 0

    amU = aL + (aU - aL)/phi
    ErmU = err_check(amU,coord_vector,slope,space_type)
    while (abs(ErmU - ErmL) > 0.0000001) and (aU - aL > 0.000000001): #While errors haven't converged and while aL and aU
        if ErmU < ErmL:
            aL = amL
            ErL = ErmL #unnecessary?
            amL = amU
            ErmL = ErmU
            amU = aL + (aU - aL)/phi
            ErmU = err_check(amU,coord_vector,slope,space_type)
            if amL > amU:
                amL,amU = amU,amL
                ErmL,ErmU = ErmU,ErmL
            #print 'amU = ' + str(amU) + '  ErmU = ' + str(ErmU)
            #print '%s %s %s %s' % (aL, amL, amU, aU)
            #print '%s %s %s %s' % (ErL, ErmL, ErmU, ErU)
        else:
            aU = amU
            ErU = ErmU #unnecessary?
            amU = amL
            ErmU = ErmL
            amL = aU - (aU - aL)/phi
            ErmL = err_check(amL,coord_vector,slope,space_type)
            if amL > amU:
                amL,amU = amU,amL
                ErmL,ErmU = ErmU,ErmL
            #print 'amL = ' + str(amL) + '  ErmL = ' + str(ErmL)
            #print '%s %s %s %s' % (aL, amL, amU, aU)
            #print '%s %s %s %s' % (ErL, ErmL, ErmU, ErU)
    
    if ErL == min(ErL,ErmL,ErmU,ErU):
        print 'returning: aL = %s, ErL = %s' % (aL,ErL)
        return aL
    elif ErmL == min(ErmL,ErmU,ErU):
        print 'returning: amL = %s, ErmL = %s' % (amL,ErmL)
        return amL
    elif ErmU == min(ErmU,ErU):
        print 'returning: amU = %s, ErmU = %s' % (amU,ErmU)
        return amU
    else: #ErU == min(Errors)
        print 'returning: aU = %s, ErU = %s' % (aU,ErU)
        return aU


def slope_find(coord_vector,space_type):
    print 'entering slope_find'
    if space_type == 'coeff':
        slope = np.zeros(coeff_terms)
        terms = coeff_terms
    elif space_type == 'valley':
        slope = np.zeros(valley_terms)
        terms = valley_terms
    else:
        raise ValueError('space_type must be "coeff" or "valley"')
    error_before = err_check(0,coord_vector,slope,space_type)
    shell = np.zeros((terms*2,terms+1))
    for i in range(terms):
        shell[2*i,i] = 1
        shell[2*i+1,i] = -1
        slope = shell[2*i,0:-1]
        shell[2*i,-1] = err_check(-1*slopeCheckAlpha,coord_vector,slope,space_type) - error_before
        slope = shell[2*i+1,0:-1]
        shell[2*i+1,-1] = err_check(-1*slopeCheckAlpha,coord_vector,slope,space_type) - error_before
    print 'shell:'
    print shell
    for i in range(terms*2):
        shell[i,:] = shell[i,:]*shell[i,-1]
    slope = sum(shell[:,0:terms])/2
    print slope

    print 'exiting slope_find'
    return slope

def slope_sign_find(coord_vector):
    plus = err_check(slopeCheckAlpha,coord_vector,slope,'coeff')
    minus = err_check(-1*slopeCheckAlpha,coord_vector,slope,'coeff')
    if plus < minus:
        return 1
    else:
        return -1

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

search_disp = 0.1
def valley_basis_find(lower_space_terms,higher_space_terms):
    global slope
    R = np.zeros((higher_space_terms,higher_space_terms))
    for i in range(higher_space_terms):
        R[i,i] = search_disp
        coeff_vector = R[i] + origin_displacement
        slope = slope_find(coeff_vector,'coeff')
        alpha = alpha_find(coeff_vector,10.0,'coeff')
        R[i] -= alpha*slope
    norms = [0.0]*higher_space_terms
    for i in range(higher_space_terms):
        norms[i] = np.linalg.norm(R[i])
    smallest_norm_index = np.argmin(norms)
    R = np.delete(R,smallest_norm_index,0)
    V = gs(R)
    return V.T

#__________________________________________________________
#coeff_vector = np.array([-3.79745533 ,-7.60036733 ,-3.80799175, -0.19730373]) #Lowest minimum found
coeff_vector = np.array([-4.0 ,-8.0 ,-4.0, 0.0])
coeff_terms = len(coeff_vector)
valley_terms = coeff_terms - 1
best_coeff_vector = coeff_vector.copy()
best_coeff_vector_error = np.inf
best_coeff_iter = 0

for prim_iter in range(100):
    print '\nbeginning primary iteration\n'
    print 'initial coeff_vector:\n'
    print coeff_vector
    print 'initial error = %s' % err_check(0,coeff_vector,np.zeros(coeff_terms),'coeff')
    
    print '\nFinding valley\n'
    slope = slope_find(coeff_vector,'coeff')
    alpha = alpha_find(coeff_vector,alpha_step,'coeff')
    coeff_vector = coeff_vector - alpha*slope
    origin_displacement = coeff_vector.copy()
    
    print '\nFinding valley basis\n'
    A = valley_basis_find(valley_terms,coeff_terms)
    valley_vector = np.zeros(valley_terms)
    
    for sec_iter in range(100):
        print '\nbeginning secondary iteration\n'
        print 'Optimizing through valley space\n'
        slope = slope_find(valley_vector,'valley')
        alpha = alpha_find(valley_vector,alpha_step,'valley')
        if alpha == 0:
            print '\nReached minimum in tangent-valley space. Searching for true valley space.\n'
            break
        valley_vector = valley_vector - alpha*slope
        print 'new valley vector:'
        print valley_vector
        print 'new error = %s' % err_check(0,valley_vector,np.zeros(valley_terms),'valley')
    
    coeff_vector = VtoC(valley_vector)
    if not np.all(coeff_vector==best_coeff_vector):
        print coeff_vector
        print best_coeff_vector
        print 'Warning!!! Losing track of best coeff_vector'
    
    if best_coeff_iter != prim_iter:
        print '\nNo best_coeff found this primary iteration. Returning best_coeff_vector'
        print 'best_coeff_vector:'
        print best_coeff_vector
        print 'Er_stat = %s' % best_coeff_vector_error
        """
        for thrd_it in range(100):
            innerValley_terms = valley_terms - 1
            B = valley_basis_find(innerValley_terms,valley_terms)
        """



        break
    
    coeff_vector = best_coeff_vector.copy()