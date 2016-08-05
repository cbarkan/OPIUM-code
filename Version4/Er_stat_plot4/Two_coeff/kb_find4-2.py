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
    f = 0.0
    for i in range(len(coeff_vector)):
        f += coeff_vector[i]*(x**i)
    return f

def fourier_cos(x,coeff_vector): #Evaluates F(x) where F is a fourier cosine series with the given coefficients
    f = 0.0
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

def err_check(alpha,coeff_vector,slope): #Runs OPIUM, checks test config errors and returns their sum
    coeff_vector_test = coeff_vector.copy()
    coeff_vector_test[jj] = coeff_vector[jj] - slope*alpha*dC
    grid_test = fit_to_func(grid,coeff_vector_test,'fourier_cos')
    setKB(grid_test)
    
    call('./opium ' + filename + ' ' + filename + '.log all rpt', shell=True)

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

def fit_to_func(grid,coeff_vector_test,form):
    if form=='fourier_cos':
        for i in range(grid_size):
            grid[i,1] = fourier_cos(x_min + (i+0.5)/resolution,coeff_vector_test)
    elif form=='poly':
        for i in range(grid_size):
            grid[i,1] = poly(x_min + (i+0.5)/resolution,coeff_vector_test)
    else:
        return 'Error: functional form not properly specified'
    return grid

phi = 1.6180339887
def alpha_find(coeff_vector,alpha_step):
    #aL / ErL == alpha_lower / Error_lower
    #amL / ErmL == alpha_middle-lower / Error_middle-lower
    #amU / ErmU == alpha_middle-upper / Error_middle-upper
    #aL / ErU == alpha_upper / Error_upper
    
    aL = slopeCheckAlpha
    amL = slopeCheckAlpha
    ErL = err_check(aL,coeff_vector,slope)
    ErmL = ErL
    print ErL
    
    #Let alpha travel
    aU = alpha_step
    ErU = err_check(aU,coeff_vector,slope)
    print ErU
    
    if ErU < ErmL:
        while ErU < ErmL:
            ErL = ErmL
            aL = amL
            ErmL = ErU
            amL = aU
            aU += alpha_step
            ErU = err_check(aU,coeff_vector,slope)
            print ErU
    else:
        print 'Alpha didnt travel, reducing alpha_step'
        count = 0
        for i in range(15):
            aU = aU/4
            alpha_step = alpha_step/4
            ErU = err_check(aU,coeff_vector,slope)
            if ErU < ErmL:
                while ErU < ErmL:
                    ErL = ErmL
                    aL = amL
                    ErmL = ErU
                    amL = aU
                    aU += alpha_step
                    ErU = err_check(aU,coeff_vector,slope)
                    print ErU
                break
            elif count >= 10:
                print 'No alpha found: returning alpha = %s' % slopeCheckAlpha
                return slopeCheckAlpha #If loop completes, no good alpha found, return 0
            else:
                count +=1

    amU = aL + (aU - aL)/phi
    ErmU = err_check(amU,coeff_vector,slope)
    while (abs(ErmU - ErmL) > 0.001) and (aU - aL > 0.0001): #While errors haven't converged and while aL and aU
        if ErmU < ErmL:
            aL = amL
            ErL = ErmL #unnecessary?
            amL = amU
            ErmL = ErmU
            amU = aL + (aU - aL)/phi
            ErmU = err_check(amU,coeff_vector,slope)
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
            ErmL = err_check(amL,coeff_vector,slope)
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

def slope_find(deltaC):
    before = err_check(0,coeff_vector,0)
    coeff_vector[jj] = coeff_vector[jj] + deltaC
    up = err_check(0,coeff_vector,0)
    up_slope = (up - before)/dC
    print 'up_slope = %s' % up_slope
    coeff_vector[jj] = coeff_vector[jj] - 2*deltaC
    down = err_check(0,coeff_vector,0)
    down_slope = (before - down)/dC
    print 'down_slope = %s' % down_slope
    if up_slope*down_slope >= 0: #If slopes have the same sign
        return (up_slope + down_slope)/2
    else:
        slope = slope_find(deltaC/2)
        return slope

#______________________________________________________________________

resolution = 50.0 #boxes per a.u.
x_min = 0.01
x_max = 1.5
coeff_vector = np.array([0.0])
terms = 8
test_range = 3.0
test_step = 0.25
dC = 0.05
slopeCheckAlpha = 0.0001
alpha_step = 5.0

grid_size = int((x_max - x_min)*resolution)
grid = np.zeros((grid_size+1,2))
grid[0,0] = np.rint(np.log(x_min/(a*(z**(-1.0/3.0))))/b + 1)
for i in np.arange(1,grid_size+1): #Set box bounds
    grid[i,0] = i*(x_max-x_min)/grid_size #UNITS: a.u.
    grid[i,0] = np.log(grid[i,0]/(a*(z**(-1.0/3.0))))/b + 1 #Converts to g.p units
    grid[i,0] = np.rint(grid[i,0]) #Ensures integer values of g.p

for jj in range(terms):
    slope = slope_find(dC)
    alpha = alpha = alpha_find(coeff_vector,alpha_step)
    coeff_vector[jj] = coeff_vector[jj] - slope*alpha*dC
    print 'old coeff_vector:'
    print coeff_vector
    coeff_vector = np.append(coeff_vector, np.array([0.0]))

x = np.linspace(x_min,x_max,grid_size+1)
w = (x_max - x_min)/(grid_size+2)
plt.bar(x,grid[:,1], width = w)
#plt.title('Atomic #: ' + str(z) +'   Er_stat=' + str(Er_stat))
plt.show()