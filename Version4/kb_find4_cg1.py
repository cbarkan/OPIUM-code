#Version 4-3: slope is an array--it can affect all coefficients
##################################
#INPUTS SECTION: Adjust the variables below as desired.
##################################
filename = 'v' #(must be string) Write this as string without a file extension (i.e. 'Fe' is correct, 'Fe.param' is incorrect)
z = 23 #Atomic number of element
tconfig_num = 5 #(must be integer) Number of test configurations, this number must match the corresponding number in the .param file
alpha_step = 5.0 #(must be float, if integer is desired, write 2.0, for example) Initial step size in the search for optimal alpha
local_orb = 's' #(must be string, either 's', 'p', 'd', or 'f') local orbital setting for augmentation function. 
a = 0.0001 #Parameter which defines grid point spacing. Do not change unless you make a corresponding change in OPIUM
b = 0.013 #Parameter which defines grid point spacing. Do not change unless you make a corresponding change in OPIUM
timeout_sec = 30
##################################
##################################

import numpy as np
from subprocess import Popen
import threading
import matplotlib.pylab as plt
fig_num = 0

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

def err_check(coeff_vector): #Runs OPIUM, checks test config errors and returns their sum
    #print 'entering err_check'
    #print coeff_vector
    #print slope
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


def err_check_slope(alpha,coeff_vector,slope): #Runs OPIUM, checks test config errors and returns their sum
    #print 'entering err_check'
    #print coeff_vector
    #print slope
    coeff_vector = coeff_vector - slope*alpha
    grid = fit_to_func(coeff_vector,'fourier_cos')
    coeff_vector = coeff_vector + slope*alpha
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

def slope_find(coeff_vector):
    print 'entering slope_find'
    print coeff_vector
    slope = np.zeros(terms)
    print slope
    error_before = err_check_slope(0,coeff_vector,slope)
    
    circle_is_necessary = False
    shell = np.zeros((terms*2,terms+1))
    for i in range(terms):
        shell[2*i,i] = 1
        shell[2*i+1,i] = -1
        slope = shell[2*i,0:terms]
        shell[2*i,-1] = err_check_slope(dC,coeff_vector,slope) - error_before
        slope = shell[2*i+1,0:terms]
        shell[2*i+1,-1] = err_check_slope(dC,coeff_vector,slope) - error_before
        #if shell[2*i,-1]*shell[2*i+1,-1] > 0:
        #    circle_is_necessary = True
    print 'shell'
    print shell

    if circle_is_necessary:
        print 'Inconsistent slopes: Circle is necessary'
        inconsistent_slope_count[0] += 1.0
        print 'Total inconsistent slopes so far: %s' % inconsistent_slope_count[0]
        slope = circle_slope_search(shell,dC)
    else:
        for i in range(terms*2):
            shell[i,:] = shell[i,-1]*shell[i,:]
        slope = sum(shell[:,0:terms])/(2*dC)
        print slope
    #Check if slope is good
    """
    slope_length = (sum(np.square(slope)))**(0.5)
    if slope_length != 0:
        error_now = err_check(slopeCheckAlpha/slope_length,coeff_vector,slope)
        print error_now
        if error_now >= error_before: #slope is bad
            slope = np.zeros(terms)
            slope_length = 0
    else: #Slope is zero length
        slope = np.zeros(terms) #Ensure slope is exactly zero for condition in script
        slope_length = 0
    """
    print 'exiting slope_find'
    return slope

def circle_slope_search(shell,dC):
    #ONLY WORKS FOR terms==2
    print 'entering circle_slope_search'
    if terms != 2:
        return np.zeros(terms)
    
    shell_points = 64 #MUST BE DIVISIBLE BY 4!!!
    big_shell = np.zeros((shell_points,terms+1))
    
    for i in range(len(big_shell)):
        big_shell[i,0] = np.cos(2*np.pi*i/shell_points)
        big_shell[i,1] = np.sin(2*np.pi*i/shell_points)
        for jj in range(len(shell)):
            if np.all(shell[jj,0:2] == big_shell[i,0:2]):
                big_shell[i,2] = shell[jj,2]
            else:
                slope = np.array([big_shell[i,0],big_shell[i,1],0])
                big_shell[i,2] = err_check(1.0,coeff_vector,slope)
    
    print big_shell
    minIndex = np.argmin(big_shell[:,2])
    slope = np.array([big_shell[minIndex,0],big_shell[minIndex,1],0])
    return slope


#______________________________________________________________________

resolution = 50.0 #boxes per a.u.
x_min = 0.01
x_max = 1.5
coeff_vector = np.array([-4.0,-8.0,-4.0])
new_terms = 1
test_range = 3.0
test_step = 0.25
dC = 0.00001
slopeCheckAlpha = 1
alpha_step = 1000.0
inconsistent_slope_count = np.array([0.0])

grid_size = int((x_max - x_min)*resolution)
grid = np.zeros((grid_size+1,2))
grid[0,0] = np.rint(np.log(x_min/(a*(z**(-1.0/3.0))))/b + 1)
for i in np.arange(1,grid_size+1): #Set box bounds
    grid[i,0] = i*(x_max-x_min)/grid_size #UNITS: a.u.
    grid[i,0] = np.log(grid[i,0]/(a*(z**(-1.0/3.0))))/b + 1 #Converts to g.p units
    grid[i,0] = np.rint(grid[i,0]) #Ensures integer values of g.p


from scipy import optimize
terms = len(coeff_vector)
for jj in range(new_terms+1):
    print 'entering jj'
    opt_coeffs = optimize.fmin_cg(err_check, coeff_vector, fprime = slope_find,epsilon = 0.00001)
    coeff_vector = np.append(opt_coeffs, np.array([0.0]))
    print coeff_vector

print 'final coeff_vector'
print coeff_vector[0:-1]
print 'final error = %s' % err_check(coeff_vector[0:-1])
print 'Total inconsistent slopes: %s' % inconsistent_slope_count[0]
x = np.linspace(x_min,x_max,grid_size+1)
w = (x_max - x_min)/(grid_size+2)
plt.bar(x,grid[:,1], width = w)
#plt.title('Atomic #: ' + str(z) +'   Er_stat=' + str(Er_stat))
plt.show()