##################################
#INPUTS SECTION: Adjust the variables below as desired.
##################################
filename = 'mn' #(must be string) Write this as string without a file extension (i.e. 'Fe' is correct, 'Fe.param' is incorrect)
initial_box_heights = [0.0,0.0] #Write initial box heights in a list with format [height1,height2,height3] where heightx is a float. This list can be any length, and it's length determines the initial number of boxes
box_splits = 1 #(must be integer) Number of times KB boxes will split in half 
max_F = 1.5 #(must be float) Outer limit of KB boxes in a.u.
local_orb = 's' #(must be string, either 's', 'p', 'd', or 'f') local orbital setting for augmentation function.
##################################
##################################

#Other input variables: Typically these should not be changed!
alpha_step = 10000.0 #(must be float, if integer is desired, write 2.0, for example) Initial step size in the search for optimal alpha
dV = 0.0000001 #(must be float) d(KB box height) for slope calculation and alpha search.
a = 0.0001 #Parameter which defines grid point spacing. Do not change unless you make a corresponding change in OPIUM
b = 0.013 #Parameter which defines grid point spacing. Do not change unless you make a corresponding change in OPIUM
slopeCheckAlpha = 1.0
loop_end_it = 50
timeout_sec = 30

import numpy as np
from subprocess import Popen
import threading
from random import uniform
import time

#Define functions:

element_list = {'H':1,'He':2,'Li':3,'Be':4,'B':5,'C':6,'N':7,'O':8,'F':9,'Ne':10,'Na':11,'Mg':12,'Al':13,'Si':14,'P':15,'S':16,'Cl':17,'Ar':18,'K':19,'Ca':20,'Sc':21,'Ti':22, 'V':23,'Cr':24,'Mn':25,'Fe':26,'Co':27,'Ni':28,'Cu':29,'Zn':30,'Ga':31,'Ge':32,'As':33,'Se':34,'Br':35,'Kr':36,'Rb':37,'Sr':38,'Y':39,'Zr':40,'Nb':41,'Mo':42,'Tc':43,'Ru':44,'Rh':45,'Pd':46,'Ag':47,'Cd':48,'In':49,'Sn':50,'Sb':51,'Te':52,'I':53,'Xe':54,'Cs':55,'Ba':56,'La':57,'Ce':58,'Pr':59,'Nd':60,'Pm':61,'Sm':62,'Eu':63,'Gd':64,'Tb':65,'Dy':66,'Ho':67,'Er':68,'Tm':69,'Yb':70,'Lu':71,'Hf':72,'Ta':73,'W':74,'Re':75,'Os':76,'Ir':77,'Pt':78,'Au':79,'Hg':80,'Tl':81,'Pb':82,'Bi':83,'Po':84,'At':85,'Rn':86,'Fr':87,'Ra':88,'Ac':89,'Th':90,'Pa':91,'U':92,'Np':93,'Pu':94,'Am':95,'Cm':96,'Bk':97,'Cf':98,'Es':99,'Fm':100,'Md':101,'No':102,'Lr':103,'Rf':104,'Db':105,'Sg':106,'Bh':107,'Hs':108,'Mt':109,'Ds':110,'Rg':111,'Cn':112}
def getParams():
    param_file = open(filename + '.param', 'r')
    param = param_file.readlines()
    param_file.close

    index = 0
    got_element = 0
    got_tconfig = 0
    num_lines = len(param)
    for jj in range(num_lines):
        if got_tconfig + got_element == 2:
            break
        elif index >= num_lines - 1:
            if got_element == 0:
                print 'ERROR: element symbol not found in .param file. Check .param file formatting.'
            if got_tconfig == 0:
                print 'ERROR: number of test configurations not found in .param file. Check .param file formatting.'
            print 'Exiting kb_find due to improper .param file format.'
            quit()
        else:
            line = param[index]
            if line.startswith('[Atom]'):
                element = param[index+1][0:2]
                if element[1] == ' ':
                    element = element[0]
                count = 1
                for i in element_list:
                    count += 1
                    if i == element:
                        z = element_list[element]
                        index+=2
                        got_element = 1
                        break
                    elif count >= len(element_list):
                        print 'Element not found: in the .param file, the first letter of the element symbol must be upper case, and the second letter lower case. Also, the kb_find dictionary of elements goes up to atomic number 112, if you are using a heavier element, add it to the dictionary element_list. If neither of these fixes the problem, check the format of the .param file.'
                        print 'Exiting kb_find'
                        quit()
            elif line.startswith('[Configs]'):
                tconfig_num = int(param[index+1][0])
                index+=2
                got_tconfig = 1
            else:
                index += 1

    return z, tconfig_num

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
    
    slope_length = (sum(np.square(slope)))**(0.5)
    if slope_length != 0:
        print 'slope_length = %s' % slope_length
        error_now = err_check(slopeCheckAlpha/slope_length,grid_test,slope)
        if error_now >= error_before: #If no good slope, then try other method of calculating slope
            print 'No slope yet: trying alternative method'
            slope = np.zeros((grid_size+1,1))
            #error_before = err_check(0,grid_test,slope)
            for i in range(grid_size): #Find slope for each bin
                grid_test[i,1] = grid_test[i,1] + dV #Change height of ith box by dV
                error_now = err_check(0,grid_test,slope) #Check results
                if error_now == float('inf'):
                    print 'got inf during slope calc'
                    continue
                else:
                    slope[i] = (error_now - error_before)/dV #Record slope
            
            error_now = err_check(slopeCheckAlpha,grid_test,slope)
            if error_now >= error_before: #If no good slope
                slope = np.zeros((grid_size+1,1))
    
    print 'slope:'
    print slope
    return slope

phi = 1.6180339887
def alpha_find(grid,alpha_step):
    #aL / ErL == alpha_lower / Error_lower
    #amL / ErmL == alpha_middle-lower / Error_middle-lower
    #amU / ErmU == alpha_middle-upper / Error_middle-upper
    #aL / ErU == alpha_upper / Error_upper
    slope_length = (sum(np.square(slope)))**(0.5)
    aL = slopeCheckAlpha/slope_length
    amL = slopeCheckAlpha/slope_length
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
        while alpha_step >= slopeCheckAlpha/slope_length:
            aU = aU/4
            alpha_step = alpha_step/4
            print 'alpha_step = %s' % alpha_step
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

#_______________________________________________________________

z, tconfig_num = getParams()

print 'filename = %s' % filename
print 'z = %s' % z
print 'box_splits = %s' % box_splits
print 'max_F = %s' % max_F
print 'tconfig_num = %s' % tconfig_num
print 'alpha_step = %s' % alpha_step
print 'local_orb = %s' % local_orb
print 'dV = %s' % dV
print 'a = %s' % a
print 'b = %s' % b
print 'slopeCheckAlpha = %s' % slopeCheckAlpha
print 'loop_end_it = %s' %loop_end_it
print 'timeout_sec = %s' %timeout_sec

#Initialize grid for KB boxes. Creates grid in grid point units
grid_size = len(initial_box_heights)
grid = np.zeros((grid_size+1,2))
grid[0,0] = np.rint(np.log(0.01/(a*(z**(-1.0/3.0))))/b + 1)
for i in np.arange(1,grid_size+1): #Set box bounds
    grid[i,0] = i*max_F/grid_size #UNITS: a.u.
    grid[i,0] = np.log(grid[i,0]/(a*(z**(-1.0/3.0))))/b + 1 #Converts to g.p units
    grid[i,0] = np.rint(grid[i,0]) #Ensures integer values of g.p

initial_box_heights = np.concatenate((np.array(initial_box_heights),[0.0]),axis = 0)
grid[:,1] = initial_box_heights
print 'initial grid:'
print grid

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
            
    #split boxes
    print 'splitting boxes \n'
    grid_size = grid_size*2
    grid_old = grid.copy()
    grid = grid_split(grid_old,grid_size)

grid_size = grid_size/2 #Sets grid size to size of grid_old to find final error
slope = np.ones((grid_size+1,1))
Er_stat = err_check(0,grid_old,slope)
print '\n' + 'final Er_stat = ' + str(Er_stat)
print 'final grid:'
print grid_old
print 'Loop_ends = %s   if this number is not 0, consider increasing slopeCheckAlpha or loop_end_it' % loop_ends
import matplotlib.pylab as plt
x = np.linspace(0,max_F,grid_size+1)
w = max_F/(grid_size+2)
plt.bar(x,grid_old[:,1], width = w)
plt.title('Atomic #: ' + str(z) +'   Er_stat=' + str(Er_stat))
plt.show()
