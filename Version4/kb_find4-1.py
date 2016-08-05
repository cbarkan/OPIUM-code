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
##################################
##################################


import numpy as np
from subprocess import call
from random import uniform
import matplotlib.pylab as plt
import math
fig_num = 0

#Define functional forms of augmentation function
def poly(x,coeff_vector): #Evaluates F(x) where F is a polynomial with the given coefficients
    f = 0
    for i in range(len(coeff_vector)):
        f += coeff_vector[i]*(x**i)
    return f

def fourier_cos(x,coeff_vector): #Evaluates F(x) where F is a fourier cosine series with the given coefficients
    f = 0
    L = x_max
    for i in range(len(coeff_vector)):
        f += coeff_vector[i]*math.cos(x*i*math.pi/L)
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

def err_check(grid): #Runs OPIUM, checks test config errors and returns Er_stat
    setKB(grid)
    call('./opium ' + filename + ' ' + filename +'.log all rpt', shell=True)
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
        Er_stat = Er[0,0] #+ Er[0,1]
    return Er_stat

def fit_to_func(grid,coeff_vector,form):
    if form=='fourier_cos':
        for i in range(grid_size):
            grid[i,1] = fourier_cos(x_min + (i+0.5)/resolution,coeff_vector)
    elif form=='poly':
        for i in range(grid_size):
            grid[i,1] = poly(x_min + (i+0.5)/resolution,coeff_vector)
    else:
        return 'Error: functional form not properly specified'
    return grid

#______________________________________________________________________

resolution = 50.0 #boxes per a.u.
x_min = 0.01
x_max = 1.5
coeff_vector = np.array([0.0])
terms = 8
test_range = 3.0
test_step = 0.25


grid_size = int((x_max - x_min)*resolution)
grid = np.zeros((grid_size+1,2))
grid[0,0] = np.rint(np.log(x_min/(a*(z**(-1.0/3.0))))/b + 1)
for i in np.arange(1,grid_size+1): #Set box bounds
    grid[i,0] = i*(x_max-x_min)/grid_size #UNITS: a.u.
    grid[i,0] = np.log(grid[i,0]/(a*(z**(-1.0/3.0))))/b + 1 #Converts to g.p units
    grid[i,0] = np.rint(grid[i,0]) #Ensures integer values of g.p


for jj in range(terms):
    test_errors = np.zeros((2*test_range/test_step,2))
    test_errors[:,0] = np.arange(-test_range,test_range,test_step)
    for i in range(len(test_errors)):
        coeff_vector[jj] = test_errors[i,0]
        grid = fit_to_func(grid,coeff_vector,'fourier_cos')
        test_errors[i,1] = err_check(grid)
    print 'test_errors:'
    print test_errors
    plt.plot(test_errors[:,0],test_errors[:,1])
    plt.show()
    bestIndex = np.argmin(test_errors[:,1])
    bestError = test_errors[bestIndex,1]
    print 'bestError = %s' % bestError
    coeff_vector[jj] = test_errors[bestIndex,0]
    print 'coeff_vector:'
    print coeff_vector
    
    coeff_vector = np.append(coeff_vector, [0])
    test_range /= 2
    test_step /= 2

x = np.linspace(x_min,x_max,grid_size+1)
w = (x_max - x_min)/(grid_size+2)
plt.bar(x,grid[:,1], width = w)
#plt.title('Atomic #: ' + str(z) +'   Er_stat=' + str(Er_stat))
plt.show()