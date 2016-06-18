"""
BEFORE USING THIS SCRIPT:
make sure the [KBdesign] block in your .param file is in the following form:

[KBdesign]
s               #This can be set however you want
15             #this must be equal to the grid_size variable below
au 0.01 ... #The inner limit of the first KB box must be set to 0.01
au 0.40 ... #Any other KB boxes must NOT have their inner limit set to 0.01. Note: all KB boxes will be overwritten by this
                   script.

HOW TO USE:
Under the "Inputs" section below, set the variables to their proper values.
put this file into a folder with opium and a .param file with the format explained above
type in command line: python kb_find.py

OUTPUT:
When the script finishes running, the .param file will contain the optimized KB box values.

"""

#Inputs
filename = 'Cu_test' #Write this as string without a file extension
grid_size = 10 #Number of KB boxes (keep as integer)
max_F = 1.5 #Outer limit of KB boxes (in a.u.) (keep as float)
tconfig_num = 7 #Number of test configurations, this number must match the corresponding number in the .param file (keep as integer)
alpha_step = 1.0 #Initial step size in the search for optimal alpha (keep as float)

#Define functions:

import numpy as np
from subprocess import call

def setKB(grid): #Sets [KBdesign] in .param file to specified grid
    #Get param file text
    param_file = open(filename + '.param', 'r')
    param = param_file.readlines()
    param_file.close

    #Open param file for rewrite
    param_file = open(filename + '.param', 'w')

    for line in param:
        if line.startswith("au 0.01"):
            index = 0
            for i in np.arange(grid_size): #This adds all the new boxes
                param_file.write("au " + str(grid[index,0]) + " " + str(grid[index+1,0]) + " " + str(grid[index,1]) + "\n")
                index += 1
        elif line.startswith("au"): #This deletes all previus boxes
            0
        else: #This copies the rest of the file as is
            param_file.write(line)
            
    param_file.close
    return


def err_check(): #Runs OPIUM, checks test config errors and returns their sum
    call('./opium ' + filename + ' ' + filename +'.log all rpt', shell=True)

    rpt_file = open(filename + '.rpt', 'r')
    rpt = rpt_file.readlines()
    Er = np.ones((tconfig_num,2))
    Er_index = 0
    for line in rpt:
        if line.startswith(" AE-NL-  t"):
            Er[Er_index,0] = float(line[22:38])
            Er[Er_index,1] = float(line[42:58])
            Er_index += 1
            
    rpt_file.close

    Er_stat = sum(Er[:,0]) #Sum of test config errors
    return Er_stat

#_____________________________________________________



#Step 1
dV = 0.2

#Step 2
grid = np.zeros((grid_size+1,2))
grid[0,0] = 0.01
for i in np.arange(1,grid_size+1):
    grid[i,0] = i*max_F/grid_size

setKB(grid)
before = err_check()
print before

#Steps 3,4,5, and 6
slope = np.ones((grid_size+1,1))

for jj in range(3):
    for i in range(grid_size): #Find slope for each bin
        grid[i,1] = grid[i,1] + dV
        setKB(grid)
        now = err_check()
        slope[i] = (now - before)/dV
        grid[i,1] = grid[i,1] - dV
        
    setKB(grid) #resets KB boxes after last iteration of loop above
        
    alpha = alpha_step
    grid_orig = grid.copy()

    #Let alpha "travel" to optimal range
    grid[:,1] = grid_orig[:,1] - alpha*dV*np.transpose(slope)
    setKB(grid)
    now = err_check()
    print now
    while now < before:
        before = now
        alpha += alpha_step
        grid[:,1] = grid_orig[:,1] - alpha*dV*np.transpose(slope)
        setKB(grid)
        now = err_check()
        print now


    #Use golden ratio method to narrow in on best alpha
    alpL = max(alpha - 2*alpha_step, 0) #set lower bound to left side of where alpha "travelled", or 0 if alpha didn't travel
    alpU = alpha #set upper bound to uppermost point of travel
    print (alpU)
    phi = 1.6180339887 #golden ratio

    gridL = grid.copy() #Saves lefthand grid
    gridU = grid.copy()
    for i in range(15): #THIS CAN BE MADE MORE EFFICIENT
        alpha = alpL + (alpU - alpL)/phi
        gridU[:,1] = grid_orig[:,1] - alpha*dV*np.transpose(slope)
        setKB(gridU)
        nowU = err_check()
        alpha = alpU - (alpU - alpL)/phi
        gridL[:,1] = grid_orig[:,1] - alpha*dV*np.transpose(slope)
        setKB(gridL)
        nowL = err_check()
        if nowU > nowL:
            alpU = alpL + (alpU - alpL)/phi
            alpha = alpL + (alpU - alpL)/phi
            before = nowL
            print(str(alpL) + ' ' + str(alpU))
            print nowL
        else:
            alpL = alpU - (alpU - alpL)/phi
            alpha = alpU - (alpU - alpL)/phi
            before = nowU
            print(str(alpL) + ' ' + str(alpU))
            print nowU

    grid[:,1] = grid_orig[:,1] - alpha*dV*np.transpose(slope)
    alpha_step = alpha_step/4
    
    print alpha_step
    print 'End of jj loop'









