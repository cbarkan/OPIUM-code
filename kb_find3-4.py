##################################
#INPUTS SECTION: Adjust the variables below as desired.
##################################
filename = 'v' #(must be string) Write this as string without a file extension (i.e. 'Fe' is correct, 'Fe.param' is incorrect)
z = 23 #Atomic number of element
initial_box_num = 4 #(must be integer) Number of KB boxes at start 
box_splits = 1 #(must be integer) Number of times KB boxes qill split in half 
max_F = 1.5 #(must be float, if integer is desired, write 2.0, for example) Outer limit of KB boxes in a.u.
tconfig_num = 5 #(must be integer) Number of test configurations, this number must match the corresponding number in the .param file
alpha_step = 3.0 #(must be float, if integer is desired, write 2.0, for example) Initial step size in the search for optimal alpha
local_orb = 's' #(must be string, either 's', 'p', 'd', or 'f') local orbital setting for augmentation function. 
dV = 0.05 #(must be float) d(KB box height) for slope calculation and alpha search.
a = 0.0001 #Parameter which defines grid point spacing. Do not change unless you make a corresponding change in OPIUM
b = 0.013 #Parameter which defines grid point spacing. Do not change unless you make a corresponding change in OPIUM
##################################
##################################

import numpy as np
from subprocess import call
from random import uniform

#Define functions:
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

def err_check(alpha,grid,slope): #Runs OPIUM, checks test config errors and returns their sum
    grid_test = grid.copy()
    grid_test[:,1] = grid[:,1] - alpha*dV*np.transpose(slope)
    setKB(grid_test)
    
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

def slope_find(grid):
    print 'Entering slope_find'
    slope = np.ones((grid_size+1,1))
    error_before = err_check(0,grid,slope)
    print error_before
    for i in range(grid_size): #Find slope for each bin
        grid[i,1] = grid[i,1] + dV #Change height of ith box by dV
        error_now = err_check(0,grid,slope) #Check results
        slope[i] = (error_now - error_before)/dV #Record slope
        grid[i,1] = grid[i,1] - dV #Restore grid to previous value
    
    
    error_now = err_check(0.05,grid,slope)
    print error_now
    if error_now > error_before: #If no good slope
        slope = np.zeros((grid_size+1,1))
        print 'found no good slope'
    
    return slope

phi = 1.6180339887
def alpha_find(grid):
    #aL / ErL == alpha_lower / Error_lower
    #amL / ErmL == alpha_middle-lower / Error_middle-lower
    #amU / ErmU == alpha_middle-upper / Error_middle-upper
    #aL / ErU == alpha_upper / Error_upper
    
    aL = 0.05
    amL = 0.05
    ErL = err_check(aL,grid,slope)
    ErmL = ErL
    print ErL
    
    #Let alpha travel
    aU = alpha_step
    ErU = err_check(aU,grid,slope)
    print ErU
    
    while ErU < ErmL:
        ErL = ErmL
        aL = amL
        ErmL = ErU
        amL = aU
        aU += alpha_step
        ErU = err_check(aU,grid,slope)
        print ErU

    amU = aL + (aU - aL)/phi
    ErmU = err_check(amU,grid,slope)
    while (abs(ErmU - ErmL) > 0.01) and (aU - aL > 0.001): #While errors haven't converged and while aL and aU
        if ErmU < ErmL:
            aL = amL
            ErL = ErmL #unnecessary?
            amL = amU
            ErmL = ErmU
            amU = aL + (aU - aL)/phi
            ErmU = err_check(amU,grid,slope)
            print 'amU = ' + str(amU) + '  ErmU = ' + str(ErmL)
        else:
            aU = amU
            ErU = ErmU #unnecessary?
            amU = amL
            ErmU = ErmL
            amL = aU - (aU - aL)/phi
            ErmL = err_check(amL,grid,slope)
            print 'amL = ' + str(amL) + '  ErmL = ' + str(ErmL)
    
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
#_______________________________________________________________

#Initialize grid for KB boxes. Creates grid in grid point units
grid_size = initial_box_num
grid = np.zeros((grid_size+1,2))
grid[0,0] = np.rint(np.log(0.01/(a*(z**(-1.0/3.0))))/b + 1)
for i in np.arange(1,grid_size+1): #Set box bounds
    grid[i,0] = i*max_F/grid_size #UNITS: a.u.
    grid[i,0] = np.log(grid[i,0]/(a*(z**(-1.0/3.0))))/b + 1 #Converts to g.p units
    grid[i,0] = np.rint(grid[i,0]) #Ensures integer values of g.p

for jj in range(box_splits+1):
    slope = np.ones((grid_size+1,1))
    for i in range(3):
        slope = slope_find(grid)
        if np.all(slope == np.zeros((grid_size+1,1))): #if slope = 0
            print 'didnt find good slope'
            break
        alpha = alpha_find(grid)
        if alpha == 0:
            print 'break: alpha=0'
            break
        else:
            grid[:,1] = grid[:,1] - alpha*dV*np.transpose(slope)
            checking_the_error = err_check(0,grid,slope)
            print checking_the_error
    print 'splitting boxes \n'
    #split boxes
    grid_size = initial_box_num * (2**(jj+1))
    grid_old = grid.copy()
    grid = np.zeros((grid_size+1,2))
    grid[0,0] = np.rint(np.log(0.01/(a*(z**(-1.0/3.0))))/b + 1)
    for i in np.arange(1,grid_size+1): #Set box bounds
        grid[i,0] = i*max_F/grid_size
        grid[i,0] = np.log(grid[i,0]/(a*(z**(-1.0/3.0))))/b + 1 #Converts to g.p units
        grid[i,0] = np.rint(grid[i,0]) #Ensures integer values of g.p
    
    for i in range(len(grid_old) - 1): #Set new box to latest box settings
        grid[2*i,1] = grid_old[i,1]
        grid[2*i+1,1] = grid_old[i,1]
        
grid_size = grid_size/2 #Sets grid size to size of grid_old to find final error
Er_stat = err_check(0,grid_old,slope)
print '\n' + 'final Er_stat = ' + str(Er_stat)
print 'final grid:'
print grid_old 
