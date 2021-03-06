"""
BEFORE USING THIS SCRIPT:
make sure the [KBdesign] block in your .param file is in the following form:

[KBdesign]
s               #This can be set to any value, but the initial value will not be saved. To choose a local orbital, adjust the local_orb variable in the input section below.
15             #This can be set to any value. The script will change and update this value as necessary will the script runs.
au 0.01     #Anything can be here, but this line must exist (i.e. the previous line cannot be the last line in the file). If this is a KB box it will be overwritten. Otherwise, the line will remain unchanged.

HOW TO USE:
Under the "Inputs" section below, set the variables to the values you desire.
Put this file into a folder with opium and a .param file with the format explained above
type into the command line: python kb_find.py

OUTPUT:
As the script runs, the step-by-step output is printed to the command line. To save this output to a file, call the script with this command: python kb_find.py > history_file.txt
Note that in the file history_file.txt, OPIUM errors do not appear in the order they occured, and thus should be ignored.
When the script finishes running, the .param file will contain the optimized KB box values.

WARNING: This script will modify the .param file and will not save the previous version. If you want to keep the original .param file, save it in a different directory or with a different name.
_____________________________________________________
New in version 3-3
-Grid point units are used instead of au. This allows boxes are spaced so that no overlap occurs

"""
##################################
#INPUTS SECTION: Adjust the variables below as desired.
##################################
filename = 'v' #(must be string) Write this as string without a file extension (i.e. 'Fe' is correct, 'Fe.param' is incorrect)
z = 23 #Atomic number of element
initial_box_num = 8 #(must be integer) Number of KB boxes at start 
box_splits = 2 #(must be integer) Number of times KB boxes qill split in half 
max_F = 1.5 #(must be float, if integer is desired, write 2.0, for example) Outer limit of KB boxes in a.u.
tconfig_num = 5 #(must be integer) Number of test configurations, this number must match the corresponding number in the .param file
alpha_step = 1.0 #(must be float, if integer is desired, write 2.0, for example) Initial step size in the search for optimal alpha
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


def err_check(): #Runs OPIUM, checks test config errors and returns their sum
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

#_____________________________________________________

#Initialize grid for KB boxes. Creates grid in grid point units
grid_size = initial_box_num
grid = np.zeros((grid_size+1,2))
grid[0,0] = np.rint(np.log(0.01/(a*(z**(-1.0/3.0))))/b + 1)
for i in np.arange(1,grid_size+1): #Set box bounds
    grid[i,0] = i*max_F/grid_size #UNITS: a.u.
    grid[i,0] = np.log(grid[i,0]/(a*(z**(-1.0/3.0))))/b + 1 #Converts to g.p units
    grid[i,0] = np.rint(grid[i,0]) #Ensures integer values of g.p

print grid
setKB(grid)
before = err_check() #Run an initial trial with no boxes: no OPIUM error should result if other parameters are chosen well

#Find random starting point that does not produce errors

before = float('inf')
while before == float('inf'):
    for i in np.arange(0,grid_size):
        grid[i,1] = uniform(2,-2)
    
    setKB(grid)
    before = err_check()


print 'Initial grid:'
print grid

for mm in range(1,box_splits+2):
    slope = np.ones((grid_size+1,1))
    for jj in range(2):
        print 'Er_stat w/ initial grid: ' + str(before) + '\n'
        for i in range(grid_size): #Find slope for each bin
            grid[i,1] = grid[i,1] + dV #Change height of ith box by dV
            setKB(grid)
            now = err_check() #Check results
            slope[i] = (now - before)/dV #Record slope
            grid[i,1] = grid[i,1] - dV #Restore grid to previous value
        print slope    
        alpha = alpha_step
        grid_orig = grid.copy()
        grid_before = grid.copy()

        #Let alpha "travel" to optimal range
        grid[:,1] = grid_orig[:,1] - alpha*dV*np.transpose(slope)
        setKB(grid)
        now = err_check()
        print 'alpha: ' + str(alpha)+ '    Er_stat: ' + str(now)
        while now < before:
            grid_before = grid.copy()
            before = now
            alpha += alpha_step
            grid[:,1] = grid_orig[:,1] - alpha*dV*np.transpose(slope)
            setKB(grid)
            now = err_check()
            print 'alpha: ' + str(alpha)+ '    Er_stat: ' + str(now)


        #Use golden section search to narrow in on best alpha
        alpL = max(alpha - 2*alpha_step, 0) #set lower bound to left side of where alpha "travelled", or 0 if alpha didn't travel
        alpU = alpha #set upper bound to uppermost point of travel
        phi = 1.6180339887 #golden ratio

        gridL = grid.copy() #Saves lefthand grid
        gridU = grid.copy()
        
        count = 0
        nowU = float('inf')
        nowL = -1*float('inf')
        while (abs(nowU - nowL) > 0.1) and (count <= 30): #THIS CAN BE MADE MORE EFFICIENT
            alpha = alpL + (alpU - alpL)/phi
            gridU[:,1] = grid_orig[:,1] - alpha*dV*np.transpose(slope)
            setKB(gridU)
            nowU = err_check()
            print '\n' + 'nowU: ' + str(nowU) +  '   alpha = ' + str(alpha)
            alpha = alpU - (alpU - alpL)/phi
            gridL[:,1] = grid_orig[:,1] - alpha*dV*np.transpose(slope)
            setKB(gridL)
            nowL = err_check()
            if nowL == float('inf'): #If gridL produced OPIUM error, ensure lower bound on alpha does not change
                nowL = -1*nowL
            print 'nowL: ' + str(nowL) +  '   alpha = ' + str(alpha)

            if (nowU < nowL) and (nowU < before):
                alpL = alpU - (alpU - alpL)/phi
                print 'new alpha bounds: ' + str(alpL) + '  ' + str(alpU)
                print gridU
            else:
                alpU = alpL + (alpU - alpL)/phi
                print 'new alpha bounds: ' + str(alpL) + '  ' + str(alpU)
                print gridL
                
            count += 1
        
        #Save results of last iteration for the return to the top of jj loop
        if (before < nowU) and (before < abs(nowL)): #No improvement from previous jj iteration
            print 'choosing before. before = ' + str(before)
            print 'End of jj loop:'
            grid = grid_before.copy()
            setKB(grid)
            Er_stat = err_check()
            print 'final grid:'
            print grid
            print 'final Er_stat'
            print Er_stat
            print '\n \n'
            break #leave loop to split boxes
        elif nowU < nowL:
            grid = gridU.copy()
            before = nowU
        else:
            grid = gridL.copy()
            before = nowL
        
        print 'End of jj loop:'
        setKB(grid)
        Er_stat = err_check()
        print 'final grid:'
        print grid
        print 'final Er_stat'
        print Er_stat
        print '\n \n'
        
    #Split grid: each box splits in half
    grid_size = initial_box_num * (2**mm)
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
    
grid_size = grid_size/2
setKB(grid_old)
Er_stat = err_check()
print '\n' + 'final Er_stat = ' + str(Er_stat)
print 'final grid:'
print grid_old