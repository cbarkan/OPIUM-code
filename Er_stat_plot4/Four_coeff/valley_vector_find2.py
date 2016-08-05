import numpy as np
import matplotlib.pylab as plt
from subprocess import Popen
import threading

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

#Define functional forms of augmentation function
def poly(x,coeff_vector): #Evaluates F(x) where F is a polynomial with the given coefficients
    f = 0.0
    for i in range(coeff_terms):
        f += coeff_vector[i]*(x**i)
    return f

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

#Define python functions
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

def err_check(alpha,valley_vector,slope): #Runs OPIUM, checks test config errors and returns their sum
    valley_vector = valley_vector - alpha*slope
    coeff_vector = VtoC(valley_vector)
    valley_vector = valley_vector + alpha*slope
    grid = fit_to_func(coeff_vector,'fourier_cos')
    setKB(grid)
    
    #Use threading to cancel OPIUM if it takes too long
    proc = Popen('./opium ' + filename + ' ' + filename +'.log all rpt', shell=True)
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
    return Er_stat

phi = 1.6180339887
def alpha_find(alpha_step):
    #aL / ErL == alpha_lower / Error_lower
    #amL / ErmL == alpha_middle-lower / Error_middle-lower
    #amU / ErmU == alpha_middle-upper / Error_middle-upper
    #aL / ErU == alpha_upper / Error_upper
    slope_length = np.linalg.norm(slope)
    if slope_length==0:
        print 'Exiting alpha_find: slope_length == 0. Returning alpha = 0'
        return 0
    aL = 0
    amL = 0
    ErL = err_check(aL,valley_vector,slope)
    ErmL = ErL
    print 'testing1:'
    print valley_vector
    print ErL
    
    #Let alpha travel
    aU = alpha_step
    ErU = err_check(aU,valley_vector,slope)
    print ErU
    
    if ErU < ErmL:
        while ErU < ErmL:
            ErL = ErmL
            aL = amL
            ErmL = ErU
            amL = aU
            aU += alpha_step
            ErU = err_check(aU,valley_vector,slope)
            print ErU
    else:
        print 'Alpha didnt travel, reducing alpha_step'
        for i in range(10000):
            aU = aU/4
            alpha_step = alpha_step/4
            print 'alpha_step = %s' % alpha_step
            ErU = err_check(aU,valley_vector,slope)
            print ErU
            if ErU < ErmL:
                print ErU
                while ErU < ErmL:
                    ErL = ErmL
                    aL = amL
                    ErmL = ErU
                    amL = aU
                    aU += alpha_step
                    ErU = err_check(aU,valley_vector,slope)
                    print ErU
                break
            elif alpha_step < 0.000000001*slopeCheckAlpha/slope_length:
                print 'No alpha found: returning alpha = %s' % 0
                return 0 #If loop completes, no good alpha found, return 0

    amU = aL + (aU - aL)/phi
    ErmU = err_check(amU,valley_vector,slope)
    while (abs(ErmU - ErmL) > 0.001) and (aU - aL > 0.0001): #While errors haven't converged and while aL and aU
        if ErmU < ErmL:
            aL = amL
            ErL = ErmL #unnecessary?
            amL = amU
            ErmL = ErmU
            amU = aL + (aU - aL)/phi
            ErmU = err_check(amU,valley_vector,slope)
            if amL > amU:
                amL,amU = amU,amL
                ErmL,ErmU = ErmU,ErmL
            print 'amU = ' + str(amU) + '  ErmU = ' + str(ErmU)
            print '%s %s %s %s' % (aL, amL, amU, aU)
            print '%s %s %s %s' % (ErL, ErmL, ErmU, ErU)
        else:
            aU = amU
            ErU = ErmU #unnecessary?
            amU = amL
            ErmU = ErmL
            amL = aU - (aU - aL)/phi
            ErmL = err_check(amL,valley_vector,slope)
            if amL > amU:
                amL,amU = amU,amL
                ErmL,ErmU = ErmU,ErmL
            print 'amL = ' + str(amL) + '  ErmL = ' + str(ErmL)
            print '%s %s %s %s' % (aL, amL, amU, aU)
            print '%s %s %s %s' % (ErL, ErmL, ErmU, ErU)
    
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

def slope_find():
    print 'entering slope_find'
    slope = np.zeros(valley_terms)
    print valley_vector
    print valley_terms
    print slope
    error_before = err_check(0,valley_vector,slope)
    shell = np.zeros((valley_terms*2,valley_terms+1))
    for i in range(valley_terms):
        shell[2*i,i] = 1
        shell[2*i+1,i] = -1
        slope = shell[2*i,0:-1]
        shell[2*i,-1] = err_check(-1*slopeCheckAlpha,valley_vector,slope) - error_before
        slope = shell[2*i+1,0:-1]
        shell[2*i+1,-1] = err_check(-1*slopeCheckAlpha,valley_vector,slope) - error_before
    print 'shell:'
    print shell
    for i in range(valley_terms*2):
        shell[i,:] = shell[i,:]*shell[i,-1]
    slope = sum(shell[:,0:valley_terms])/2
    print slope

    print 'exiting slope_find'
    return slope

#_________________________________________________________________________

plotname = 'zoom_n4-n8-n4-0'
mesh = np.load(plotname + '_Eig1.npy') + np.load(plotname + '_Norm1.npy')
c0_index = np.load(plotname + '_c0_index.npy')
c1_index = np.load(plotname + '_c1_index.npy')
c2_index = np.load(plotname + '_c2_index.npy')
c3_index = np.load(plotname + '_c3_index.npy')
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


#_________________________________________________________________________________________

###Inputs###
filename = 'v' #(must be string) Write this as string without a file extension (i.e. 'Fe' is correct, 'Fe.param' is incorrect)
z = 23 #Atomic number of element
tconfig_num = 5 #(must be integer) Number of test configurations, this number must match the corresponding number in the .param file
alpha_step = 5.0 #(must be float, if integer is desired, write 2.0, for example) Initial step size in the search for optimal alpha
local_orb = 's' #(must be string, either 's', 'p', 'd', or 'f') local orbital setting for augmentation function. 
a = 0.0001 #Parameter which defines grid point spacing. Do not change unless you make a corresponding change in OPIUM
b = 0.013 #Parameter which defines grid point spacing. Do not change unless you make a corresponding change in OPIUM
timeout_sec = 30

alpha_step = 1000.0
slopeCheckAlpha = 0.00002
resolution = 50.0 #boxes per a.u.
x_min = 0.01
x_max = 1.5
valley_vector = np.array([ 0.38729833,  0.31515762 , 0.22511258])
new_terms = 0
############

grid_size = int((x_max - x_min)*resolution)
grid = np.zeros((grid_size+1,2))
grid[0,0] = np.rint(np.log(x_min/(a*(z**(-1.0/3.0))))/b + 1)
for i in np.arange(1,grid_size+1): #Set box bounds
    grid[i,0] = i*(x_max-x_min)/grid_size #UNITS: a.u.
    grid[i,0] = np.log(grid[i,0]/(a*(z**(-1.0/3.0))))/b + 1 #Converts to g.p units
    grid[i,0] = np.rint(grid[i,0]) #Ensures integer values of g.p

coeff_terms = 4
valley_terms = 3

print 'initital valley_vector:'
print valley_vector
print 'initial error = %s' % err_check(0,valley_vector,np.zeros(3))
for i in range(3): 
    slope = slope_find()
    alpha = alpha_find(alpha_step)
    valley_vector = valley_vector - alpha*slope
    print 'valley_vector:'
    print valley_vector
    print 'coeff_vector:'
    print VtoC(valley_vector)
    print 'error = %s' % err_check(0,valley_vector,np.zeros(3))








"""
x = np.array([ 0.38729833,  0.31515762 , 0.22511258])
coeff_vect = VtoC(x)
print coeff_vect
a,b,c,d = CtoInd(coeff_vect)
print a,b,c,d
print mesh[a,b,c,d]
"""
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