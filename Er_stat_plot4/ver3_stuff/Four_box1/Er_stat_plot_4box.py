##################################
#INPUTS SECTION: Adjust the variables below as desired.
##################################
filename = 'v' #(must be string) Write this as string without a file extension (i.e. 'Fe' is correct, 'Fe.param' is incorrect)
z = 23 #Atomic number of element
initial_box_num = 2 #(must be integer) Number of KB boxes at start 
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
timeout_sec = 30
##################################
##################################


import numpy as np
from subprocess import Popen
import threading
import matplotlib.pylab as plt
from matplotlib.colors import LogNorm
import time

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

def err_check(grid): #Runs OPIUM, checks test config errors and returns their sum
    #THIS IS NOT THE STANDARD err_check FUNCTION!!!!!
    grid_test = grid.copy()
    setKB(grid_test)
    
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
        Er_eig = float('inf') #Returns inf so that current KB box is not chosen
        Er_norm = float('inf')
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
        #Er_stat = Er[0,0] + Er[0,1]
        Er_eig = Er[0,0]
        Er_norm = Er[0,1]
    return Er_eig, Er_norm


grid_size = 4
grid = np.zeros((grid_size+1,2))
grid[0,0] = np.rint(np.log(0.01/(a*(z**(-1.0/3.0))))/b + 1)
for i in np.arange(1,grid_size+1): #Set box bounds
    grid[i,0] = i*max_F/grid_size #UNITS: a.u.
    grid[i,0] = np.log(grid[i,0]/(a*(z**(-1.0/3.0))))/b + 1 #Converts to g.p units
    grid[i,0] = np.rint(grid[i,0]) #Ensures integer values of g.p


box0_lower_bound = -10
box0_upper_bound = 0
box1_lower_bound = -10
box1_upper_bound = 0
box2_lower_bound = -5
box2_upper_bound = 5
box3_lower_bound = -5
box3_upper_bound = 5
mesh_step = 1.0
Eig_mesh = np.zeros(((box0_upper_bound - box0_lower_bound)/mesh_step,(box1_upper_bound - box1_lower_bound)/mesh_step,(box2_upper_bound - box2_lower_bound)/mesh_step,(box3_upper_bound - box3_lower_bound)/mesh_step,))
Norm_mesh = np.zeros(((box0_upper_bound - box0_lower_bound)/mesh_step,(box1_upper_bound - box1_lower_bound)/mesh_step,(box2_upper_bound - box2_lower_bound)/mesh_step,(box3_upper_bound - box3_lower_bound)/mesh_step,))
box0_index = np.zeros((box0_upper_bound - box0_lower_bound)/mesh_step)
box1_index = np.zeros((box1_upper_bound - box1_lower_bound)/mesh_step)
box2_index = np.zeros((box2_upper_bound - box2_lower_bound)/mesh_step)
box3_index = np.zeros((box3_upper_bound - box3_lower_bound)/mesh_step)

start_time = time.time()

for i0 in range(int((box0_upper_bound - box0_lower_bound)/mesh_step)):
    grid[0,1] = i0*mesh_step + box0_lower_bound
    box0_index[i0] = i0*mesh_step + box0_lower_bound
    for i1 in range(int((box1_upper_bound - box1_lower_bound)/mesh_step)):
        grid[1,1] = i1*mesh_step + box1_lower_bound
        box1_index[i1] = i1*mesh_step + box1_lower_bound
        for i2 in range(int((box2_upper_bound - box2_lower_bound)/mesh_step)):
            grid[2,1] = i2*mesh_step + box2_lower_bound
            box2_index[i2] = i2*mesh_step + box2_lower_bound
            for i3 in range(int((box3_upper_bound - box3_lower_bound)/mesh_step)):
                grid[3,1] = i3*mesh_step + box3_lower_bound
                box3_index[i3] = i3*mesh_step + box3_lower_bound
                print grid
                #if np.all(grid[:,1] == np.array([-5.0,-2.0,3.0,-4.0,0])):
                   # continue
                Er_eig, Er_norm = err_check(grid)
                print 'Er_eig = %s   Er_norm = %s' % (Er_eig, Er_norm)
                Eig_mesh[i0,i1,i2,i3] = Er_eig
                Norm_mesh[i0,i1,i2,i3] = Er_norm
elapsed_time = time.time() - start_time
print 'Elapsed time: %s' % elapsed_time

Total_mesh = Eig_mesh + Norm_mesh


run_name = 'n10-0-n5-5'

"""
plt.figure(0)
plt.imshow(Eig_mesh, interpolation='none', origin = 'lower',extent=[Right_lower_bound,Right_upper_bound,Left_lower_bound,Left_upper_bound], norm=LogNorm(vmin=1, vmax=5000))
cb = plt.colorbar(orientation='horizontal')
cb.set_label('Eigenvalue error (mRy)')
plt.ylabel('Inner Box height (0.01 to 0.75 au)')
plt.xlabel('Outer Box height (0.75 to 1.5 au)')
plt.savefig(str(run_name) + '_eig.png')

plt.figure(1)
plt.imshow(Norm_mesh, interpolation='none', origin = 'lower',extent=[Right_lower_bound,Right_upper_bound,Left_lower_bound,Left_upper_bound], norm=LogNorm(vmin=1, vmax=5000))
cb = plt.colorbar(orientation='horizontal')
cb.set_label('Tail norm error (mRy)')
plt.ylabel('Inner Box height (0.01 to 0.75 au)')
plt.xlabel('Outer Box height (0.75 to 1.5 au)')
plt.savefig(str(run_name) + '_norm.png')

plt.figure(2)
plt.imshow(Total_mesh, interpolation='none', origin = 'lower',extent=[Right_lower_bound,Right_upper_bound,Left_lower_bound,Left_upper_bound], norm=LogNorm(vmin=1, vmax=5000))
cb = plt.colorbar(orientation='horizontal')
cb.set_label('Eigenvalue error + Norm error(mRy)')
plt.ylabel('Inner Box height (0.01 to 0.75 au)')
plt.xlabel('Outer Box height (0.75 to 1.5 au)')
plt.savefig(str(run_name) + '_total.png')
"""
np.save(str(run_name) + '_Eig.npy', Eig_mesh)
np.save(str(run_name) + '_Norm.npy', Norm_mesh)
np.save(str(run_name) + '_Total.npy',Total_mesh)
np.save(str(run_name) + '_box0_index.npy', box0_index)
np.save(str(run_name) + '_box1_index.npy', box1_index)
np.save(str(run_name) + '_box2_index.npy', box2_index)
np.save(str(run_name) + '_box3_index.npy', box3_index)
