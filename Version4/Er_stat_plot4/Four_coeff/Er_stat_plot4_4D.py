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
resolution = 50.0 #boxes per a.u.
x_min = 0.01
x_max = 1.5
terms = 8
test_range = 3.0
test_step = 0.25
dC = 0.05
slopeCheckAlpha = 0.0001
alpha_step = 5.0
##################################
##################################
timeout_sec = 30

import numpy as np
from subprocess import Popen
import threading
import matplotlib.pylab as plt
import math
import time
from matplotlib.colors import LogNorm
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

def err_check(coeff_vector): #Runs OPIUM, checks test config errors and returns their sum
    coeff_vector_test = coeff_vector.copy()
    grid_test = fit_to_func(grid,coeff_vector_test,'fourier_cos')
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
        Er_eig1,Er_norm1,Er_eig2,Er_norm2,Er_eig3,Er_norm3,Er_eig4,Er_norm4,Er_eig5,Er_norm5 = float('inf')*np.ones(10)
        DD = np.ones(tconfig_num*(tconfig_num+1)/2)*float('inf')
    else: #No error occured
        Er = np.ones((tconfig_num,2))
        Er_index = 0
        DD = np.ones(tconfig_num*(tconfig_num+1)/2)
        DD_index = 0
        for line in rpt:
            if line.startswith(" AE-NL-  t"):
                Er[Er_index,0] = float(line[22:38])
                Er[Er_index,1] = float(line[42:58])
                Er_index += 1
            elif line.startswith(" AE-NL-   ") and (not line.startswith(" AE-NL-   i")):
                DD[DD_index] = abs(float(line[22:33]))
                DD_index += 1
        
        rpt_file.write('#kb_find read-receipt')
        rpt_file.close

        #Er_stat = sum(Er[:,0]) #Sum of test config eigenvalue errors
        #Er_stat = sum(sum(DD)) #Sum of energy differences
        #Er_stat = sum(sum(Er)) #Sum of eigenvalue errors and norm errors
        #Er_stat = sum(sum(Er))/tconfig_num + sum(sum(DD))/(tconfig_num*(tconfig_num-1)/2) #Sum of avg eigenvalue error, avg norm error, and avg energy difference
        #Er_stat = Er[0,0] #+ Er[0,1]
        Er_eig1 = Er[0,0]
        Er_norm1 = Er[0,1]
        Er_eig2 = Er[1,0]
        Er_norm2 = Er[1,1]
        Er_eig3 = Er[2,0]
        Er_norm3 = Er[2,1]
        Er_eig4 = Er[3,0]
        Er_norm4 = Er[3,1]
        Er_eig5 = Er[4,0]
        Er_norm5 = Er[4,1]
    return Er_eig1, Er_norm1, Er_eig2, Er_norm2, Er_eig3, Er_norm3, Er_eig4, Er_norm4, Er_eig5, Er_norm5, DD

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
#___________________________________________________________________________

grid_size = int((x_max - x_min)*resolution)
grid = np.zeros((grid_size+1,2))
grid[0,0] = np.rint(np.log(x_min/(a*(z**(-1.0/3.0))))/b + 1)
for i in np.arange(1,grid_size+1): #Set box bounds
    grid[i,0] = i*(x_max-x_min)/grid_size #UNITS: a.u.
    grid[i,0] = np.log(grid[i,0]/(a*(z**(-1.0/3.0))))/b + 1 #Converts to g.p units
    grid[i,0] = np.rint(grid[i,0]) #Ensures integer values of g.p

coeff_vector = np.array([0.0,0.0,0.0,0.0])

c0_lower_bound = -4.5
c0_upper_bound = -3.5
c1_lower_bound = -8.5
c1_upper_bound = -7.5
c2_lower_bound = -4.5
c2_upper_bound = -3.5
c3_lower_bound = -0.5
c3_upper_bound = 0.5
mesh_step = 0.1

"""
c0_lower_bound = -1
c0_upper_bound = 1
c1_lower_bound = -1
c1_upper_bound = 1
c2_lower_bound = -1
c2_upper_bound = 1
c3_lower_bound = -1
c3_upper_bound = 1
mesh_step = 1.0
"""

coeff0_trials = int((c0_upper_bound - c0_lower_bound)/mesh_step)
coeff1_trials = int((c1_upper_bound - c1_lower_bound)/mesh_step)
coeff2_trials = int((c2_upper_bound - c2_lower_bound)/mesh_step)
coeff3_trials = int((c3_upper_bound - c3_lower_bound)/mesh_step)
num_trials = coeff0_trials*coeff1_trials*coeff2_trials*coeff3_trials

Eig1_mesh = np.zeros((coeff0_trials,coeff1_trials,coeff2_trials,coeff3_trials))
Norm1_mesh = Eig1_mesh.copy()
Eig2_mesh = Eig1_mesh.copy()
Norm2_mesh = Norm1_mesh.copy()
Eig3_mesh = Eig1_mesh.copy()
Norm3_mesh = Norm1_mesh.copy()
Eig4_mesh = Eig1_mesh.copy()
Norm4_mesh = Norm1_mesh.copy()
Eig5_mesh = Eig1_mesh.copy()
Norm5_mesh = Norm1_mesh.copy()
DD_mesh = np.ones((coeff0_trials,coeff1_trials,coeff2_trials,coeff3_trials,(tconfig_num*(tconfig_num+1)/2)))

c0_index = np.zeros(coeff0_trials)
c1_index = np.zeros(coeff1_trials)
c2_index = np.zeros(coeff2_trials)
c3_index = np.zeros(coeff3_trials)

start_time = time.time()
current_trial = 0
print num_trials
for c0 in range(coeff0_trials):
    coeff_vector[0] = c0*mesh_step + c0_lower_bound
    c0_index[c0] = c0*mesh_step + c0_lower_bound
    for c1 in range(coeff1_trials):
        coeff_vector[1] = c1*mesh_step + c1_lower_bound
        c1_index[c1] = c1*mesh_step + c1_lower_bound
        for c2 in range(coeff2_trials):
            coeff_vector[2] = c2*mesh_step + c2_lower_bound
            c2_index[c2] = c2*mesh_step + c2_lower_bound
            for c3 in range(coeff3_trials):
                coeff_vector[3] = c3*mesh_step + c3_lower_bound
                c3_index[c3] = c3*mesh_step + c3_lower_bound
                print coeff_vector
                Er_eig1, Er_norm1, Er_eig2, Er_norm2, Er_eig3, Er_norm3, Er_eig4, Er_norm4, Er_eig5, Er_norm5, DD = err_check(coeff_vector)
                
                Eig1_mesh[c0,c1,c2,c3] = Er_eig1
                Norm1_mesh[c0,c1,c2,c3] = Er_norm1
                Eig2_mesh[c0,c1,c2,c3] = Er_eig2
                Norm2_mesh[c0,c1,c2] = Er_norm2
                Eig3_mesh[c0,c1,c2,c3] = Er_eig3
                Norm3_mesh[c0,c1,c2,c3] = Er_norm3
                Eig4_mesh[c0,c1,c2,c3] = Er_eig4
                Norm4_mesh[c0,c1,c2,c3] = Er_norm4
                Eig5_mesh[c0,c1,c2,c3] = Er_eig5
                Norm5_mesh[c0,c1,c2,c3] = Er_norm5
                DD_mesh[c0,c1,c2,c3:] = DD
                current_trial += 1
            time_so_far = (time.time() - start_time)/60
            ETA = (time_so_far/current_trial)*(num_trials - current_trial)
            print 'Time so far: %s' % time_so_far
            print 'Estimated minutes until completion: %s' % ETA
elapsed_time = time.time() - start_time
print 'Elapsed time: %s' % elapsed_time
time_per_trial = elapsed_time/num_trials
print 'Time per trial: %s' % time_per_trial

run_name = 'zoom_n4-n8-n4-0'
"""
plt.figure(0)
plt.imshow(Eig_mesh, interpolation='none', origin = 'lower',extent=[c1_lower_bound,c1_upper_bound,c0_lower_bound,c0_upper_bound], norm=LogNorm(vmin=1, vmax=5000))
cb = plt.colorbar(orientation='horizontal')
cb.set_label('Eigenvalue error (mRy)')
plt.ylabel('Inner Box height (0.01 to 0.75 au)')
plt.xlabel('Outer Box height (0.75 to 1.5 au)')
plt.savefig(str(run_name) + '_eig.png')

plt.figure(1)
plt.imshow(Norm_mesh, interpolation='none', origin = 'lower',extent=[c1_lower_bound,c1_upper_bound,c0_lower_bound,c0_upper_bound], norm=LogNorm(vmin=1, vmax=5000))
cb = plt.colorbar(orientation='horizontal')
cb.set_label('Tail norm error (mRy)')
plt.ylabel('Inner Box height (0.01 to 0.75 au)')
plt.xlabel('Outer Box height (0.75 to 1.5 au)')
plt.savefig(str(run_name) + '_norm.png')

plt.figure(2)
plt.imshow(Total_mesh, interpolation='none', origin = 'lower',extent=[c1_lower_bound,c1_upper_bound,c0_lower_bound,c0_upper_bound], norm=LogNorm(vmin=1, vmax=5000))
cb = plt.colorbar(orientation='horizontal')
cb.set_label('Eigenvalue error + Norm error(mRy)')
plt.ylabel('Inner Box height (0.01 to 0.75 au)')
plt.xlabel('Outer Box height (0.75 to 1.5 au)')
plt.savefig(str(run_name) + '_total.png')
"""
np.save(str(run_name) + '_Eig1.npy', Eig1_mesh)
np.save(str(run_name) + '_Norm1.npy', Norm1_mesh)
np.save(str(run_name) + '_Eig2.npy', Eig2_mesh)
np.save(str(run_name) + '_Norm2.npy', Norm2_mesh)
np.save(str(run_name) + '_Eig3.npy', Eig3_mesh)
np.save(str(run_name) + '_Norm3.npy', Norm3_mesh)
np.save(str(run_name) + '_Eig4.npy', Eig4_mesh)
np.save(str(run_name) + '_Norm4.npy', Norm4_mesh)
np.save(str(run_name) + '_Eig5.npy', Eig5_mesh)
np.save(str(run_name) + '_Norm5.npy', Norm5_mesh)
np.save(str(run_name) + '_DD_mesh.npy', DD_mesh)
np.save(str(run_name) + '_c0_index.npy', c0_index)
np.save(str(run_name) + '_c1_index.npy', c1_index)
np.save(str(run_name) + '_c2_index.npy', c2_index)
np.save(str(run_name) + '_c3_index.npy', c3_index)