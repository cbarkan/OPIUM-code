##################################
#INPUTS SECTION: Adjust the variables below as desired.
##################################
filename = 'v' #(must be string) Write this as string without a file extension (i.e. 'Fe' is correct, 'Fe.param' is incorrect)
z = 23 #Atomic number of element
initial_box_num = 4 #(must be integer) Number of KB boxes at start 
box_splits = 1 #(must be integer) Number of times KB boxes qill split in half 
max_F = 1.5 #(must be float, if integer is desired, write 2.0, for example) Outer limit of KB boxes in a.u.
tconfig_num = 5 #(must be integer) Number of test configurations, this number must match the corresponding number in the .param file
alpha_step = 5.0 #(must be float, if integer is desired, write 2.0, for example) Initial step size in the search for optimal alpha
local_orb = 's' #(must be string, either 's', 'p', 'd', or 'f') local orbital setting for augmentation function. 
dV = 0.05 #(must be float) d(KB box height) for slope calculation and alpha search.
a = 0.0001 #Parameter which defines grid point spacing. Do not change unless you make a corresponding change in OPIUM
b = 0.013 #Parameter which defines grid point spacing. Do not change unless you make a corresponding change in OPIUM
slopeCheckAlpha = 0.0001
loop_end_it = 50
weight = 0.5 #How to weight average for box_smooth
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
    grid_test = grid.copy()
    slope = np.zeros((grid_size+1,1))
    error_before = err_check(0,grid_test,slope)
    for i in range(grid_size): #Find slope for each bin
        grid_test[i,1] = grid_test[i,1] + dV #Change height of ith box by dV
        error_now = err_check(0,grid_test,slope) #Check results
        if error_now == float('inf'):
            grid_test[i,1] = grid_test[i,1] - dV #Restore grid to previous value
            print 'got inf during slope calc'
        else:
            slope[i] = (error_now - error_before)/dV #Record slope
            grid_test[i,1] = grid_test[i,1] - dV #Restore grid to previous value
    
    error_now = err_check(slopeCheckAlpha,grid_test,slope)
    if error_now > error_before: #If no good slope, then try other method of calculating slope
        print 'No slope yet: trying alternative method'
        slope = np.zeros((grid_size+1,1))
        error_before = err_check(0,grid_test,slope)
        for i in range(grid_size): #Find slope for each bin
            grid_test[i,1] = grid_test[i,1] + dV #Change height of ith box by dV
            error_now = err_check(0,grid_test,slope) #Check results
            if error_now == float('inf'):
                print 'got inf during slope calc'
                continue
            else:
                slope[i] = (error_now - error_before)/dV #Record slope
        
        error_now = err_check(slopeCheckAlpha,grid_test,slope)
        if error_now > error_before: #If no good slope
            slope = np.zeros((grid_size+1,1))
    
    print 'slope:'
    print slope
    return slope

def slope_find_shift(grid,boxes_per_block,shift):
    grid_test = grid.copy()
    block_num = (grid_size - shift) / boxes_per_block #rounds down because all variables are int
    slope = np.zeros((grid_size+1,1))
    error_before = err_check(0,grid_test,slope)
        
    #Left end
    grid_test[0:shift,1] = grid_test[0:shift,1] + dV
    error_now = err_check(0,grid_test,slope) #Check results
    grid_test[0:shift,1] = grid_test[0:shift,1] - dV #Restore grid to previous value
    if error_now == float('inf'):
        print 'got inf during slope calc'
    else:
        slope[0:shift] = (error_now - error_before)/dV #Record slope
        
    #Middle blocks
    for i in range(block_num):
        left_box = shift+i*boxes_per_block
        right_box = left_box + boxes_per_block
        grid_test[left_box:right_box,1] = grid_test[left_box:right_box,1] + dV
        error_now = err_check(0,grid_test,slope) #Check results
        grid_test[left_box:right_box,1] = grid_test[left_box:right_box,1] - dV #Restore grid to previous value
        if error_now == float('inf'):
            print 'got inf during slope calc'
        else:
            slope[left_box:right_box] = (error_now - error_before)/dV #Record slope

    #Right end
    grid_test[right_box+1:,1] = grid_test[right_box+1:,1] + dV
    error_now = err_check(0,grid_test,slope) #Check results
    grid_test[right_box+1:,1] = grid_test[right_box+1:,1] - dV
    if error_now == float('inf'):
        print 'got inf during slope calc'
    else:
        slope[right_box+1:] = (error_now - error_before)/dV #Record slope
    
    print 'shifted slope:'
    print slope
    
    error_now = err_check(slopeCheckAlpha,grid_test,slope)
    if error_now > error_before: #If no good slope
        slope = np.zeros((grid_size+1,1))
    
    print 'shifted slope:'
    print slope
    return slope

phi = 1.6180339887
def alpha_find(grid,alpha_step):
    #aL / ErL == alpha_lower / Error_lower
    #amL / ErmL == alpha_middle-lower / Error_middle-lower
    #amU / ErmU == alpha_middle-upper / Error_middle-upper
    #aL / ErU == alpha_upper / Error_upper
    
    aL = slopeCheckAlpha
    amL = slopeCheckAlpha
    ErL = err_check(aL,grid,slope)
    ErmL = ErL
    print ErL
    
    #Let alpha travel
    aU = alpha_step
    ErU = err_check(aU,grid,slope)
    print ErU
    
    if ErU < ErmL:
        while ErU < ErmL:
            ErL = ErmL
            aL = amL
            ErmL = ErU
            amL = aU
            aU += alpha_step
            ErU = err_check(aU,grid,slope)
            print ErU
    else:
        print 'Alpha didnt travel, reducing alpha_step'
        count = 0
        for i in range(15):
            aU = aU/4
            alpha_step = alpha_step/4
            ErU = err_check(aU,grid,slope)
            if ErU < ErmL:
                while ErU < ErmL:
                    ErL = ErmL
                    aL = amL
                    ErmL = ErU
                    amL = aU
                    aU += alpha_step
                    ErU = err_check(aU,grid,slope)
                    print ErU
                break
            elif count >= 10:
                print 'No alpha found: returning alpha = %s' % slopeCheckAlpha
                return slopeCheckAlpha #If loop completes, no good alpha found, return 0
            else:
                count +=1

    amU = aL + (aU - aL)/phi
    ErmU = err_check(amU,grid,slope)
    while (abs(ErmU - ErmL) > 0.001) and (aU - aL > 0.0001): #While errors haven't converged and while aL and aU
        if ErmU < ErmL:
            aL = amL
            ErL = ErmL #unnecessary?
            amL = amU
            ErmL = ErmU
            amU = aL + (aU - aL)/phi
            ErmU = err_check(amU,grid,slope)
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
            ErmL = err_check(amL,grid,slope)
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
    
def grid_split(grid_old,grid_size):
    grid = np.zeros((grid_size+1,2))
    grid[0,0] = np.rint(np.log(0.01/(a*(z**(-1.0/3.0))))/b + 1)
    for i in np.arange(1,grid_size+1): #Set box bounds
        grid[i,0] = i*max_F/grid_size
        grid[i,0] = np.log(grid[i,0]/(a*(z**(-1.0/3.0))))/b + 1 #Converts to g.p units
        grid[i,0] = np.rint(grid[i,0]) #Ensures integer values of g.p
    
    for i in range(len(grid_old) - 1): #Set new box to latest box settings
        grid[2*i,1] = grid_old[i,1]
        grid[2*i+1,1] = grid_old[i,1]
        
    return grid

def box_smooth(weight):
    slope = np.ones((grid_size+1,1))
    err_notsmooth = err_check(0,grid,slope)
    print 'not smooth grid:'
    print grid
    print 'err_notsmooth = ' + str(err_notsmooth)
    
    grid_smooth = grid.copy()
    grid_smooth[0,1] = grid_old[0,1]
    grid_smooth[-2,1] = grid_old[-2,1]
    for i in range(1,len(grid_old) - 1): #Set new box to latest box settings
        grid_smooth[2*i-1,1] = grid_old[i-1,1]*(weight) + grid_old[i,1]*(1 - weight)
        grid_smooth[2*i,1] = grid_old[i-1,1]*(weight) + grid_old[i,1]*(1 -weight)
    err_smooth = err_check(0,grid_smooth,slope)
    print 'smooth grid:'
    print grid_smooth
    print 'err_smooth = ' + str(err_smooth)
    if err_smooth == float('inf'):
        grid_smooth = grid
        print 'Smooth grid gave OPIUM error'

    return grid_smooth

def adjust_individual_box(grid,test_range,test_step,fig_num):
    grid_test = grid.copy()
    bestErrors = np.ones((grid_size,2))
    for box in range(grid_size):
        test_errors = np.zeros((2*test_range/test_step,2))
        for mm in range(int(2*test_range/test_step)):
            height = -1*test_range + mm*test_step + grid[box,1]
            grid_test[box,1] = height
            test_errors[mm,0] = height
            test_errors[mm,1] = err_check(0,grid_test,slope)
        grid_test[box,1] = grid[box,1]
        print 'test_errors for box %s of %s' % (box,grid_size-1)
        print test_errors

        bestIndex = np.argmin(test_errors[:,1])
        bestErrors[box,:] = test_errors[bestIndex,:]

        plt.figure(fig_num)
        fig_num += 1
        plt.plot(test_errors[:,0],test_errors[:,1])
        title = 'box_%s_of_%s' % (box,grid_size-1)
        plt.title('Error vs Height,' + title)
        plt.savefig(title + '.png')

    print 'bestErrors; the bestHeight and its error for each box:'
    print bestErrors
    bestIndex = np.argmin(bestErrors[:,1])
    print 'chose box %s' % bestIndex
    grid_test[bestIndex,1] = bestErrors[bestIndex,0]
    print 'grid:'
    print grid_test
    return grid_test, fig_num

#_______________________________________________________________

print 'filename = %s' % filename
print 'z = %s' % z
print 'initial_box_num = %s' % initial_box_num
print 'box_splits = %s' % box_splits
print 'max_F = %s' % max_F
print 'tconfig_num = %s' % tconfig_num
print 'alpha_step = %s' % alpha_step
print 'local_orb = %s' % local_orb
print 'dV = %s' % dV
print 'a = %s' % a
print 'b = %s' % b
print 'slopeCheckAlpha = %s' % slopeCheckAlpha

#Initialize grid for KB boxes. Creates grid in grid point units
grid_size = initial_box_num
grid = np.zeros((grid_size+1,2))
grid[0,0] = np.rint(np.log(0.01/(a*(z**(-1.0/3.0))))/b + 1)
for i in np.arange(1,grid_size+1): #Set box bounds
    grid[i,0] = i*max_F/grid_size #UNITS: a.u.
    grid[i,0] = np.log(grid[i,0]/(a*(z**(-1.0/3.0))))/b + 1 #Converts to g.p units
    grid[i,0] = np.rint(grid[i,0]) #Ensures integer values of g.p

grid[0,1] = -6.0
grid[1,1] = -4.0
grid[2,1] = 3.0
grid[3,1] = 0.0

loop_ends = 0 #Number of times for-loop to recalculate slope ended without reaching minimum
for jj in range(box_splits+1):
    slope = np.ones((grid_size+1,1))
    count = 0
    for i in range(loop_end_it):
        slope = slope_find(grid)
        if np.all(slope == np.zeros((grid_size+1,1))): #if slope = 0
            print 'No good slope found'
            break
        alpha = alpha_find(grid,alpha_step)
        grid[:,1] = grid[:,1] - alpha*dV*np.transpose(slope)
        print grid
        count +=1
        if count >= loop_end_it - 1:
            loop_ends += 1
            
    #split and smooth boxes
    print 'splitting boxes \n'
    grid_size = initial_box_num * (2**(jj+1))
    grid_old = grid.copy()
    grid = grid_split(grid_old,grid_size)
    #grid = box_smooth(weight)

#grid_size = grid_size/2 #Sets grid size to size of grid_old to find final error
grid_size = grid_size/2
slope = np.ones((grid_size+1,1))
Er_stat = err_check(0,grid_old,slope)
print '\n' + 'final Er_stat = ' + str(Er_stat)
print 'final grid:'
print grid_old
print 'Loop_ends = %s   if this number is not 0, increase loop_end_it' % loop_ends
import matplotlib.pylab as plt
x = np.linspace(0,max_F,grid_size+1)
w = max_F/(grid_size+2)
plt.bar(x,grid_old[:,1], width = w)
plt.title('Atomic #: ' + str(z) +'   Er_stat=' + str(Er_stat))
plt.show()

#___________________________________________________________________
#Additional code:
"""

#This loop runs through optimization with slope_shift
    for i in range(loop_end_it):
        slope = slope_find_shift(grid,2,1)
        if np.all(slope == np.zeros((grid_size+1,1))): #if slope = 0
            print 'No good slope found'
            break
        alpha = alpha_find(grid)
        grid[:,1] = grid[:,1] - alpha*dV*np.transpose(slope)
        print grid



#This code runs through a version2 search
    print 'entering ver2 search'
    reruns = 2
    fineness = 2.0
    slope = np.ones((grid_size+1,1))
    now = err_check(0,grid,slope)
    for mm in np.arange(fineness + 0.1):
        for run in range(reruns):
            for i in np.arange(grid_size):
                grid[i,1] += 0.2/(2**mm)
                before = now
                now = err_check(0,grid,slope)
                if now >= before:
                    grid[i,1] -= 0.4/(2**mm)
                    now = err_check(0,grid,slope)
                    if now >= before:
                        grid[i,1] += 0.2/(2**mm)
                        now = before
                        print 'found no change for i=%s run=%s mm=%s' % (i,run,mm)
                    else:
                        print now
                        print grid
                else:
                    print now
                    print grid
    print 'exiting ver2 search'


"""
