import numpy as np
import matplotlib.pyplot as plt
import math

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


filename = 'v'
local_orb = 's'
a = 0.0001
b = 0.013
z = 23
resolution = 30.0 #boxes per a.u.
x_min = 0.01
x_max = 3.01
coeff_vector = np.array([1,.333333,0.2,1.0/7.0])

def poly(x,coeff_vector): #Evaluates F(x) where F is a polynomial with the given coefficients
    f = 0
    for i in coeff_vector:
        f += coeff_vector[i]*(x**i)
    return f

def fourier_cos(x,coeff_vector):
    f = 0
    L = x_max
    for i in coeff_vector:
        f += coeff_vector[i]*math.cos(x*i*math.pi/L)
    return f

grid_size = int((x_max - x_min)*resolution)
grid = np.zeros((grid_size+1,2))

grid[0,0] = np.rint(np.log(x_min/(a*(z**(-1.0/3.0))))/b + 1)
for i in np.arange(1,grid_size+1): #Set box bounds
    grid[i,0] = i*(x_max-x_min)/grid_size #UNITS: a.u.
    grid[i,0] = np.log(grid[i,0]/(a*(z**(-1.0/3.0))))/b + 1 #Converts to g.p units
    grid[i,0] = np.rint(grid[i,0]) #Ensures integer values of g.p


for i in range(grid_size):
    grid[i,1] = fourier_cos(x_min + (i+0.5)/res,coeff_vector)
    
x = np.linspace(x_min,x_max,grid_size+1)
w = (x_max - x_min)/(grid_size+2)
plt.bar(x,grid[:,1], width = w)
#plt.title('Atomic #: ' + str(z) +'   Er_stat=' + str(Er_stat))
plt.show()