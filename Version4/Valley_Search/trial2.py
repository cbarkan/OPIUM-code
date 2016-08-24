import numpy as np
from subprocess import Popen
import threading

###Inputs###
filename = 'v' #(must be string) Write this as string without a file extension (i.e. 'Fe' is correct, 'Fe.param' is incorrect)
z = 23 #Atomic number of element
tconfig_num = 5 #(must be integer) Number of test configurations, this number must match the corresponding number in the .param file
alpha_step = 5.0 #(must be float, if integer is desired, write 2.0, for example) Initial step size in the search for optimal alpha
local_orb = 's' #(must be string, either 's', 'p', 'd', or 'f') local orbital setting for augmentation function. 
a = 0.0001 #Parameter which defines grid point spacing. Do not change unless you make a corresponding change in OPIUM
b = 0.013 #Parameter which defines grid point spacing. Do not change unless you make a corresponding change in OPIUM
timeout_sec = 30

alpha_step = 100000.0
slopeCheckAlpha = 0.0000005
resolution = 50.0 #boxes per a.u.
x_min = 0.01
x_max = 1.5
############

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

def fourier_cos(x,coeff_vector): #Evaluates F(x) where F is a fourier cosine series with the given coefficients
	f = 0.0
	L = x_max
	for i in range(coeff_terms):
		f += coeff_vector[i]*np.cos(x*i*np.pi/L)
	return f

def fit_to_func(coeff_vector_test,form):
	if form=='fourier_cos':
		for i in range(grid_size):
			grid[i,1] = fourier_cos(x_min + (i+0.5)/resolution,coeff_vector_test)
	elif form=='poly':
		for i in range(grid_size):
			grid[i,1] = poly(x_min + (i+0.5)/resolution,coeff_vector_test)
	else:
		raise ValueError('Error: functional form not properly specified')
	return grid

def VtoC(x,level):
	for i in range(level,0,-1):
		x = np.dot(COB[i-1],x) + origin_disps[i-1]
	return x

def err_check(alpha,coord_vector,slope,level): #Runs OPIUM, checks test config errors and returns their sum
	test_coeff_vector = VtoC(coord_vector - alpha*slope,level)
	grid = fit_to_func(test_coeff_vector,'fourier_cos')
	setKB(grid)

	#Use threading to cancel OPIUM if it takes too long
	proc = Popen('./opium ' + filename + ' ' + filename + '.log all rpt', shell=True)
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
			if line.startswith("  !!ERROR!!"):
				#Ghosts in psp
				return float('inf')
			if line.startswith(" AE-NL-  t"):
				Er[Er_index,0] = float(line[22:38])
				Er[Er_index,1] = float(line[42:58])
				Er_index += 1
			elif line.startswith(" AE-NL-   ") and (not line.startswith(" AE-NL-   i")):
				DD[DD_index,0] = abs(float(line[22:33]))
				DD_index += 1
		
		rpt_file.write('#kb_find read-receipt')
		rpt_file.close

		Er_stat = Er[0,0] + Er[0,1]

		global best_coeff_vector
		global best_coeff_vector_error
		global best_cycleNum
		if Er_stat < best_coeff_vector_error:
			best_coeff_vector_error = Er_stat
			best_coeff_vector = test_coeff_vector.copy()
			best_cycleNum = cycleNum
			print 'New best coeff_vector!!'
			
	return Er_stat

phi = 1.6180339887
def alpha_find(slope,coord_vector,alpha_step,level):
	#aL / ErL == alpha_lower / Error_lower
	#amL / ErmL == alpha_middle-lower / Error_middle-lower
	#amU / ErmU == alpha_middle-upper / Error_middle-upper
	#aL / ErU == alpha_upper / Error_upper

	#Check slope. This needs more thought...
	print 'Entering alpha_find'

	slope_length = np.linalg.norm(slope)
	if slope_length==0:
		print 'Exiting alpha_find: slope_length == 0. Returning alpha = 0'
		return 0.0, 0.0
	alpha_test = 2.0
	Er0 = err_check(0,coord_vector,slope,level)
	Er1 = err_check(alpha_test,coord_vector,slope,level)
	if Er1 > Er0:
		alpha_test *= -1.0
		Er1 = err_check(alpha_test,coord_vector,slope,level)
		if Er1 < Er0:
			print 'Reversing slope direction'
			slope *= -1.0
			alpha_test *= -1.0 #change alpha_test back to original
		else:
			alpha_test /= 4
			Er2 = err_check(alpha_test,coord_vector,slope,level)
			while alpha_test > 0.001*slopeCheckAlpha/slope_length:
				alpha_test /= 4
				Er1 = err_check(alpha_test,coord_vector,slope,level)
				if Er1 < Er0:
					break
			else:
				print 'No good alpha found, returning 0'
				return 0.0,0.0

	aL = alpha_test
	amL = alpha_test
	ErL = Er1
	ErmL = ErL
	
	#Let alpha travel
	aU = alpha_step + aL
	ErU = err_check(aU,coord_vector,slope,level)
	print ErU
	
	if ErU < ErmL:
		while ErU < ErmL:
			ErL = ErmL
			aL = amL
			ErmL = ErU
			amL = aU
			aU *= 2.0
			ErU = err_check(aU,coord_vector,slope,level)
			print ErU
	else:
		print 'Alpha didnt travel, reducing alpha_step'
		while aU > alpha_test:
			aU = aU/8
			print 'alpha_step = %s' % aU
			ErU = err_check(aU,coord_vector,slope,level)
			print ErU
			if ErU < ErmL:
				print ErU
				while ErU < ErmL:
					ErL = ErmL
					aL = amL
					ErmL = ErU
					amL = aU
					aU *= 2
					ErU = err_check(aU,coord_vector,slope,level)
					print ErU
				break
		else:
			return alpha_test, Er0 - Er1 #If loop completes, no good alpha found, return 0

	amU = aL + (aU - aL)/phi
	ErmU = err_check(amU,coord_vector,slope,level)
	while abs(ErmU - ErmL) > 0.00000001: #While errors haven't converged and while aL and aU
		if ErmU < ErmL:
			aL = amL
			ErL = ErmL #unnecessary?
			amL = amU
			ErmL = ErmU
			amU = aL + (aU - aL)/phi
			ErmU = err_check(amU,coord_vector,slope,level)
			if amL > amU:
				amL,amU = amU,amL
				ErmL,ErmU = ErmU,ErmL
			#print 'amU = ' + str(amU) + '  ErmU = ' + str(ErmU)
			#print '%s %s %s %s' % (aL, amL, amU, aU)
			#print '%s %s %s %s' % (ErL, ErmL, ErmU, ErU)
		else:
			aU = amU
			ErU = ErmU #unnecessary?
			amU = amL
			ErmU = ErmL
			amL = aU - (aU - aL)/phi
			ErmL = err_check(amL,coord_vector,slope,level)
			if amL > amU:
				amL,amU = amU,amL
				ErmL,ErmU = ErmU,ErmL
			#print 'amL = ' + str(amL) + '  ErmL = ' + str(ErmL)
			#print '%s %s %s %s' % (aL, amL, amU, aU)
			#print '%s %s %s %s' % (ErL, ErmL, ErmU, ErU)
	
	if ErL == min(ErL,ErmL,ErmU,ErU):
		error_improvement = Er0 - ErL
		print 'returning: aL = %s, ErL = %s, Er_imp = %s' % (aL,ErL,error_improvement)
		return aL, error_improvement
	elif ErmL == min(ErmL,ErmU,ErU):
		error_improvement = Er0 - ErmL
		print 'returning: aL = %s, ErL = %s, Er_imp = %s' % (aL,ErL,error_improvement)
		return amL, error_improvement
	elif ErmU == min(ErmU,ErU):
		error_improvement = Er0 - ErmU
		print 'returning: aL = %s, ErL = %s, Er_imp = %s' % (aL,ErL,error_improvement)
		return amU, error_improvement
	else: #ErU == min(Errors)
		error_improvement = Er0 - ErU
		print 'returning: aL = %s, ErL = %s, Er_imp = %s' % (aL,ErL,error_improvement)
		return aU, error_improvement

def slope_find(coord_vector,level):
	print 'entering slope_find'
	terms = coeff_terms - level
	slope = np.zeros(terms)
	error_before = err_check(0,coord_vector,slope,level)
	shell = np.zeros((terms*2,terms+1))
	for i in range(terms):
		shell[2*i,i] = 1
		shell[2*i+1,i] = -1
		slope = shell[2*i,0:-1]
		shell[2*i,-1] = err_check(-1*slopeCheckAlpha,coord_vector,slope,level) - error_before
		slope = shell[2*i+1,0:-1]
		shell[2*i+1,-1] = err_check(-1*slopeCheckAlpha,coord_vector,slope,level) - error_before
	print 'shell:'
	print shell
	for i in range(terms*2):
		shell[i,:] = shell[i,:]*shell[i,-1]
	slope = sum(shell[:,0:terms])/2
	print 'slope:'
	print slope

	print 'exiting slope_find'
	return slope

def gs(X):
	#Returns matrix Y whose rows are an orthonormal basis of span{rows of X}
	Y = X.copy()
	Y[0] = Y[0]/np.linalg.norm(Y[0])
	for ii in range(1,Y.shape[0]):
		proj = 0
		for i in range(ii):
			proj += Y[i] * np.dot(Y[ii],Y[i])/np.dot(Y[i],Y[i])
		Y[ii] = Y[ii] - proj
		Y[ii] = Y[ii]/np.linalg.norm(Y[ii])
	return Y

def valley_basis_find(level):
	#level should be dim of higher space
	global slope
	search_disp = 0.1
	terms = coeff_terms - level
	R = np.zeros((terms,terms))
	for i in range(terms):
		R[i,i] = search_disp
		temp_vector = R[i] + origin_disps[level]
		slope = slope_find(temp_vector,level)
		alpha = alpha_find(slope,temp_vector,10000.0,level)[0]
		R[i] -= alpha*slope
	norms = [0.0]*terms
	for i in range(terms):
		norms[i] = np.linalg.norm(R[i])
	smallest_norm_index = np.argmin(norms)
	R = np.delete(R,smallest_norm_index,0)
	V = gs(R)
	return V.T

def cg(coord_vector,level):
	print 'Entering cg'
	#Conjugate gradient optimization algorithm from PDF on my email
	terms = coeff_terms - level
	k = 0
	r = -1*slope_find(coord_vector,level)
	d = r
	delta_new = np.dot(r,r)
	err_orig = err_check(0,coord_vector,d,level)
	while True:
		alpha, Er_imp = alpha_find(-1*d,coord_vector,alpha_step,level)
		if Er_imp < 0.0000001:
			coord_vector = coord_vector + alpha*d
			total_improvement = err_orig - err_check(0,coord_vector,d,level)
			print 'exiting cg. total_imp = %s' % total_improvement
			return coord_vector, total_improvement
		coord_vector += alpha*d
		r = -1*slope_find(coord_vector,level)
		delta_old = delta_new
		delta_new = np.dot(r,r)
		B = delta_new/delta_old
		d = r + B*d
		k += 1
		if k==terms or np.dot(r,d)<=0:
			d = r
			k = 0

###Initialization__________________________________
grid_size = int((x_max - x_min)*resolution)
grid = np.zeros((grid_size+1,2))
grid[0,0] = np.rint(np.log(x_min/(a*(z**(-1.0/3.0))))/b + 1)
for i in np.arange(1,grid_size+1): #Set box bounds
	grid[i,0] = i*(x_max-x_min)/grid_size #UNITS: a.u.
	grid[i,0] = np.log(grid[i,0]/(a*(z**(-1.0/3.0))))/b + 1 #Converts to g.p units
	grid[i,0] = np.rint(grid[i,0]) #Ensures integer values of g.p

coeff_vector = np.array([-4.0 ,-8.0 ,-4.0, 0.0])
coeff_terms = len(coeff_vector)
totalLevels = 4 #1 means no valley find. Must be >= 1 and <= coeff_terms
best_coeff_vector = coeff_vector.copy()
best_coeff_vector_error = np.inf
best_cycleNum = 0

COB = [np.zeros((coeff_terms,coeff_terms-1))]
for i in range(1,totalLevels-1):
	COB += [np.zeros((coeff_terms-i,coeff_terms-i-1))]

origin_disps = [np.zeros(coeff_terms)]
for i in range(1,totalLevels-1):
	origin_disps += [np.zeros(coeff_terms-i)]
###________________________________________________

cycleNum = 0
while True:
	level = 0
	coord_vector = best_coeff_vector.copy()
	coord_vector, Er_imp = cg(coord_vector,level)
	cycleNum += 1
	for level in range(1, totalLevels):
		print 'Entering higher level. Level=%s' % level
		print 'cycleNum=%s' % cycleNum
		origin_disps[level-1] = coord_vector.copy()
		COB[level-1] = valley_basis_find(level-1)
		coord_vector = np.zeros(coeff_terms - level)
		coord_vector, Er_imp = cg(coord_vector,level)
		if (Er_imp > 0.000001) or (best_cycleNum == cycleNum): #The best_cycleNum condition accounts for low errors found in slope calculations. Is the other condition necessary?
			break
		else:
			continue
	else:
		break

print best_coeff_vector
print best_coeff_vector_error