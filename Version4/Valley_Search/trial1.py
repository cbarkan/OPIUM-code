def cg(coord_vector,levelsDeep):
	#Conjugate gradient optimization algorithm from PDF on my email
	terms = coeff_terms - levelsDeep
	k = 0
	r = -1*slope_find()
	d = r
	delta_new = np.dot(r,r)
	err_orig = err_check(0,coord_vector,d,levelsDeep)
	while True:
		alpha, Er_imp = alpha_find(-1*d,coeff_vector,alpha_step,levelsDeep)
		if Er_imp < 0.0000001:
			coeff_vector = coeff_vector + alpha*d
			total_improvement = err_orig - err_check(0,coord_vector,d,levelsDeep)
			return coord_vector, total_improvement
		coord_vector += alpha*d
		r = -1*slope_find()
		delta_old = delta_new
		delta_new = np.dot(r,r)
		B = delta_new/delta_old
		d = r + B*d
		k += 1
		if k==terms or np.dot(r,d)<=0:
			d = r
			k = 0