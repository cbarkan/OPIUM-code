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
As the script runs, it will list the latest KBdesign settings. The last setting is automatically put into the .param file. The results can be viewed in the .rpt file.

"""
#Inputs
filename = 'Fe_test' #Write this as string without a file extension
grid_size = 10.0 #Number of KB boxes
max_F = 1.5 #Outer limit of KB boxes (in a.u.)
tconfig_num = 4 #Number of test configurations
fineness = 2.0 #Box height will be calculated to +- 2/(2^fineness)
reruns = 4

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
            Er[Er_index,0] = float(line[22:36])
            Er[Er_index,1] = float(line[42:56])
            Er_index += 1
            
    rpt_file.close

    Er_stat = sum(Er[:,0]) #Sum of test config errors
    return Er_stat



#Create grid for KB dimensions and insert it into [KBdesign] block
grid = np.zeros((grid_size+1,2))
grid[0,0] = 0.01
for i in np.arange(1,grid_size+1):
    grid[i,0] = i*max_F/grid_size

setKB(grid)

#Save grid to history
with open('box_history.txt','w') as f_handle: #Here, all contents of box_history are overwritten
    np.savetxt(f_handle,grid,delimiter=' ',footer='\n')


#Run opium and check errors
now = err_check()

#Adjust grid and check for improvements
for jj in np.arange(fineness + 0.1):
    for run in range(reruns):
        for i in np.arange(grid_size):
            grid[i,1] += 2/(2**jj)
            setKB(grid)
            before = now
            now = err_check()
            if now >= before:
                grid[i,1] -= 4/(2**jj)
                setKB(grid)
                now = err_check()
                if now >= before:
                    grid[i,1] += 2/(2**jj)
                    now = before
                else:
                    print now
                    print grid
                    with open('box_history.txt','a') as f_handle: #Here, new grid is appended to box_history
                        np.savetxt(f_handle,grid,delimiter=' ',footer='\n')
                        f_handle.write(str(jj) + ' ' + str(run) + ' ' + str(i) + '\n')
            else:
                print now
                print grid
                with open('box_history.txt','a') as f_handle: #Here, new grid is appended to box_history
                    np.savetxt(f_handle,grid,delimiter=' ',footer='\n')
                    f_handle.write(str(jj) + ' ' + str(run) + ' ' + str(i) + '\n')
