import numpy as np
from subprocess import Popen
import threading
import matplotlib.pylab as plt
from mpl_toolkits.mplot3d import Axes3D

#Define functional forms of augmentation function
def poly(x,coeff_vector): #Evaluates F(x) where F is a polynomial with the given coefficients
    f = 0.0
    for i in range(terms):
        f += coeff_vector[i]*(x**i)
    return f

def fourier_cos(x,coeff_vector): #Evaluates F(x) where F is a fourier cosine series with the given coefficients
    f = 0.0
    L = x_max
    for i in range(terms):
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

def err_check_orig(alpha,coeff_vector,d): #Runs OPIUM, checks test config errors and returns their sum
    coeff_vector = coeff_vector + d*alpha
    grid = fit_to_func(coeff_vector,'fourier_cos')
    coeff_vector = coeff_vector - d*alpha
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

test_A = np.array([[ 3.53738667, -1.83163503, -1.04345406],[-1.83163503 , 7.11373881 , 0.1887602 ],[-1.04345406,  0.1887602 ,  1.79958654]])
test_B = np.array([4,-1,3.1])
def err_check(alpha,coeff_vector,d):
    coeff_vector = coeff_vector + d*alpha
    test_error = np.dot(test_A,coeff_vector)
    test_error = 0.5*np.dot(coeff_vector,test_error)
    test_error -= np.dot(test_B,coeff_vector)
    coeff_vector = coeff_vector - d*alpha
    return test_error

phi = 1.6180339887
def alpha_find(coeff_vector,alpha_step):
    #aL / ErL == alpha_lower / Error_lower
    #amL / ErmL == alpha_middle-lower / Error_middle-lower
    #amU / ErmU == alpha_middle-upper / Error_middle-upper
    #aL / ErU == alpha_upper / Error_upper
    d_length = np.linalg.norm(d)
    aL = 0
    amL = 0
    ErL = err_check(aL,coeff_vector,d)
    ErmL = ErL
    print 'testing1:'
    print coeff_vector
    print ErL
    
    #Let alpha travel
    aU = alpha_step
    ErU = err_check(aU,coeff_vector,d)
    print ErU
    
    if ErU < ErmL:
        while ErU < ErmL:
            ErL = ErmL
            aL = amL
            ErmL = ErU
            amL = aU
            aU += alpha_step
            ErU = err_check(aU,coeff_vector,d)
            print ErU
    else:
        print 'Alpha didnt travel, reducing alpha_step'
        for i in range(100):
            aU = aU/4
            alpha_step = alpha_step/4
            print 'alpha_step = %s' % alpha_step
            ErU = err_check(aU,coeff_vector,d)
            if ErU < ErmL:
                print ErU
                while ErU < ErmL:
                    ErL = ErmL
                    aL = amL
                    ErmL = ErU
                    amL = aU
                    aU += alpha_step
                    ErU = err_check(aU,coeff_vector,d)
                    print ErU
                break
            elif alpha_step < slopeCheckAlpha/d_length:
                print 'No alpha found: returning alpha = %s' % 0
                return 0 #If loop completes, no good alpha found, return 0

    amU = aL + (aU - aL)/phi
    ErmU = err_check(amU,coeff_vector,d)
    while (abs(ErmU - ErmL) > 0.001) and (aU - aL > 0.0001): #While errors haven't converged and while aL and aU
        if ErmU < ErmL:
            aL = amL
            ErL = ErmL #unnecessary?
            amL = amU
            ErmL = ErmU
            amU = aL + (aU - aL)/phi
            ErmU = err_check(amU,coeff_vector,d)
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
            ErmL = err_check(amL,coeff_vector,d)
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
    slope = np.zeros(terms)
    print coeff_vector
    print terms
    print slope
    error_before = err_check(0,coeff_vector,slope)
    shell = np.zeros((terms*2,terms+1))
    for i in range(terms):
        shell[2*i,i] = 1
        shell[2*i+1,i] = -1
        slope = shell[2*i,0:terms]
        shell[2*i,-1] = err_check(slopeCheckAlpha,coeff_vector,slope) - error_before
        slope = shell[2*i+1,0:terms]
        shell[2*i+1,-1] = err_check(slopeCheckAlpha,coeff_vector,slope) - error_before
    print 'shell'
    print shell
    for i in range(terms*2):
        shell[i,:] = shell[i,:]*shell[i,-1]
    slope = sum(shell[:,0:terms])/2
    print slope
    print 'exiting slope_find'
    return slope

def path_plot():
    global coeff_hist
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    coeff_hist = np.vstack((coeff_hist,coeff_vector.copy()))
    ax.plot(coeff_hist[:,0],coeff_hist[:,1],coeff_hist[:,2],marker='o')
    plt.show()

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.plot(coeff_hist[-3:,0],coeff_hist[-3:,1],coeff_hist[-3:,2],marker='o')
    plt.show()

    np.save('coeff_hist.npy',coeff_hist)
#________________________________________________________________________
###Inputs###
filename = 'v' #(must be string) Write this as string without a file extension (i.e. 'Fe' is correct, 'Fe.param' is incorrect)
z = 23 #Atomic number of element
tconfig_num = 5 #(must be integer) Number of test configurations, this number must match the corresponding number in the .param file
alpha_step = 5.0 #(must be float, if integer is desired, write 2.0, for example) Initial step size in the search for optimal alpha
local_orb = 's' #(must be string, either 's', 'p', 'd', or 'f') local orbital setting for augmentation function. 
a = 0.0001 #Parameter which defines grid point spacing. Do not change unless you make a corresponding change in OPIUM
b = 0.013 #Parameter which defines grid point spacing. Do not change unless you make a corresponding change in OPIUM
timeout_sec = 30

terms = 2
alpha_step = 10000.0
slopeCheckAlpha = 0.000002
resolution = 50.0 #boxes per a.u.
x_min = 0.01
x_max = 1.5
coeff_vector = np.array([-3.5,-7.5,-3.5])
new_terms = 0
############

grid_size = int((x_max - x_min)*resolution)
grid = np.zeros((grid_size+1,2))
grid[0,0] = np.rint(np.log(x_min/(a*(z**(-1.0/3.0))))/b + 1)
for i in np.arange(1,grid_size+1): #Set box bounds
    grid[i,0] = i*(x_max-x_min)/grid_size #UNITS: a.u.
    grid[i,0] = np.log(grid[i,0]/(a*(z**(-1.0/3.0))))/b + 1 #Converts to g.p units
    grid[i,0] = np.rint(grid[i,0]) #Ensures integer values of g.p

terms = len(coeff_vector)

coeff_hist = coeff_vector.copy()


for ii in range(new_terms+1):
    k = 0
    r = -1*slope_find()
    d = r
    delta_new = np.dot(r,r)

    best_Coeffs = coeff_vector
    print 'testing2:'
    print coeff_vector
    best_Err = err_check(0,coeff_vector,np.zeros(terms))
    print 'initial error = %s' % best_Err
    for i in range(3):
        print 'starting i iteration. i = %s' % i
        print 'testing3:'
        print coeff_vector
        alpha = alpha_find(coeff_vector,alpha_step)
        coeff_vector = coeff_vector + alpha*d
        #path_plot()
        current_Err = err_check(0,coeff_vector,np.zeros(terms))
        print 'coeff_vector:'
        print coeff_vector
        print 'Error = %s' % current_Err
        if current_Err >= best_Err:
            coeff_vector = best_Coeffs
            print '\nNo improvement with new coeff_vector\n'
            break
        else:
            best_Err = current_Err
            best_Coeffs = coeff_vector
        r = -1*slope_find()
        if np.all(r==np.zeros(terms)):
            print '0 slope'
            break
        delta_old = delta_new
        delta_new = np.dot(r,r)
        B = delta_new/delta_old
        d = r + B*d
        k += 1
        if k==terms or np.dot(r,d)<=0:
            d = r
            k = 0
    
    terms += 1
    coeff_vector = np.append(coeff_vector, np.array([0.0]))

print 'true solution:'
tru_sol = np.linalg.solve(test_A,test_B)
print tru_sol
print 'true error'
print err_check(0,tru_sol,np.zeros(3))