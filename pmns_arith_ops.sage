####################################################################################################
###### IMPORTANT: parameters are supposed generated correctly so that operations are feasible ######
####################################################################################################

#~ returns a representation of 'op/phi'
def pmns_red_int(op, ri_mat, neg_iri_mat, iM_dom, phi):
	
	n = ri_mat.dimensions()[0]
	t = list(op)
	op = vector(ZZ, t+[0]*(n-len(t)))
	
	q = vector(iM_dom, op) * neg_iri_mat
	t = vector(ZZ, q) * ri_mat
	r = op + t
	return (r/phi)


#~ returns a representation of 'op/phi'
def pmns_trans_red_int(op, ri_mat, neg_iri_mat, iM_dom, phi, trans_v):
	
	n = len(trans_v)
	t = list(op)
	
	op = vector(ZZ, t+[0]*(n-len(t)))
	
	return pmns_red_int(op+trans_v, ri_mat, neg_iri_mat, iM_dom, phi)


#~ returns a representation of 'op/phi'
# ~ important: assumes phi is a power of two
def pmns_comp_red_int(op, ri_mat, neg_iri_mat, iM_dom, phi, phi_log2):
	
	n = ri_mat.dimensions()[0]
	t = list(op)
	op = vector(ZZ, t+[0]*(n-len(t)))
	
	q = vector(iM_dom, op) * neg_iri_mat
	q = vector(ZZ, q)
	i = phi_log2 - 1
	u = vector(ZZ,[(v>>i) for v in q])
	t = q * ri_mat
	w = u * ri_mat
	r = (op + t)/phi
	return (r - w)


########################################################################
########################################################################

#returns 0 if and only if A and B represent the same element.
def equality_check(A, B, ri_mat, neg_iri_mat, iM_dom, phi, trans_v):
	C = A - B 
	red = pmns_trans_red_int(C, ri_mat, neg_iri_mat, iM_dom, phi, trans_v)
	# ~ red = pmns_red_int(C, ri_mat, neg_iri_mat, iM_dom, phi)
	# ~ print(red)
	return (min(red) == max(red) == 0) 
	
########################################################################
########################################################################

#~ returns a representation of '(op1*op2)/phi'
def pmns_mont_mult(op1, op2, ri_mat, neg_iri_mat, iM_dom, phi, ext_pol, R):
	op11 = R(list(op1))
	op22 = R(list(op2))
	c = (op11 * op22)% ext_pol
	c = vector(ZZ, c)
	return pmns_red_int(c, ri_mat, neg_iri_mat, iM_dom, phi)


#~ returns a representation of '(op1*op2)/phi'
# ~ for a result in D_1 if ops are in the PMNS
def pmns_trans_mont_mult(op1, op2, ri_mat, neg_iri_mat, iM_dom, phi, ext_pol, trans_v, R):
	op11 = R(list(op1))
	op22 = R(list(op2))
	c = (op11 * op22)% ext_pol
	c = vector(ZZ, c)
	return pmns_trans_red_int(c, ri_mat, neg_iri_mat, iM_dom, phi, trans_v)


#~ returns a representation of 'op'
def exact_red_coeffs(op, ri_mat, neg_iri_mat, iM_dom, phi, ext_pol, phi2_rep, trans_v, R):
	tmp = pmns_red_int(op, ri_mat, neg_iri_mat, iM_dom, phi)
	return pmns_trans_mont_mult(tmp, phi2_rep, ri_mat, neg_iri_mat, iM_dom, phi, ext_pol, trans_v, R)

########################################################################

# ~ note: 'op' must be a PMNS element
def compute_rep_in_H(op, n, omega_poly, ri_mat, neg_iri_mat, iM_dom, phi, ext_pol, trans_v, R):
	c = pmns_trans_mont_mult(op, omega_poly, ri_mat, neg_iri_mat, iM_dom, phi, ext_pol, trans_v, R)
	for i in range(n):
		c = pmns_red_int(c, ri_mat, neg_iri_mat, iM_dom, phi)
	return c
	
# ~ note: 'op' must be a PMNS element
def compute_rep_in_HH(op, n, delta_poly, ri_mat, neg_iri_mat, iM_dom, phi, ext_pol, trans_v, phi_log2, R):
	c = pmns_trans_mont_mult(op, delta_poly, ri_mat, neg_iri_mat, iM_dom, phi, ext_pol, trans_v, R)
	for i in range(n):
		c = pmns_red_int(c, ri_mat, neg_iri_mat, iM_dom, phi)
	return pmns_comp_red_int(c, ri_mat, neg_iri_mat, iM_dom, phi, phi_log2)

########################################################################

# ~ Assumes: 'conv_polys[i]' is a rep of 'convbase^i * phi^2' or 'convbase^i * phi' according to 'in_mont_domain' value
def conv_to_pmns__to_D1(el, n, convbase_log2, phi, conv_polys, ri_mat, neg_iri_mat, iM_dom, Zp, trans_v):
	conv_mask = (1<<convbase_log2) - 1
	el = Integer(Zp(el))
	rep = 0
	for i in range(n):
		rep += (el & conv_mask) * conv_polys[i]
		el >>= convbase_log2
	return pmns_trans_red_int(rep, ri_mat, neg_iri_mat, iM_dom, phi, trans_v)


def exact_conv_to_pmns__to_D1(val, n, p, phi, ri_mat, neg_iri_mat, iM_dom, Zp, trans_v, in_mont_domain=False):
	nb_iter = n
	if in_mont_domain:
		phi_pow = phi.powermod(nb_iter + 1, p)
	else:
		phi_pow = phi.powermod(nb_iter, p)
	t = [Integer(Zp(val * phi_pow))] + [0]*(n-1)
	rep = vector(ZZ, t)
	for i in range(nb_iter-1):
		rep = pmns_red_int(rep, ri_mat, neg_iri_mat, iM_dom, phi)
	rep = pmns_trans_red_int(rep, ri_mat, neg_iri_mat, iM_dom, phi, trans_v)
	return rep


# ~ note: uses Horner approach
def conv_from_pmns(el_rep, n, gmm, phi, Zq, in_mont_domain):
	el_rep_coeffs = list(el_rep)
	el_rep_coeffs += [0]*(n - len(el_rep_coeffs))
	val = el_rep_coeffs[-1]
	for i in range(2, n+1):
		val = ((val*gmm) + el_rep_coeffs[-i])
	if in_mont_domain:
		return Integer(Zq(val/phi))
	else:
		return Integer(Zq(val))

########################################################################

def compute_conv_polys__in_D1(p, n, convbase_log2, phi, ri_mat, neg_iri_mat, iM_dom, Zp, trans_v, ext_pol, R, in_mont_domain):
	convbase = 1 << convbase_log2
	rp = exact_conv_to_pmns__to_D1(convbase*phi, n, p, phi, ri_mat, neg_iri_mat, iM_dom, Zp, trans_v,  False)
	tp = exact_conv_to_pmns__to_D1(phi, n, p, phi, ri_mat, neg_iri_mat, iM_dom, Zp, trans_v, in_mont_domain)
	conv_polys = [tp]
	for i in range(1,n):
		tp = pmns_trans_mont_mult(tp, rp, ri_mat, neg_iri_mat, iM_dom, phi, ext_pol, trans_v, R)
		conv_polys.append(tp)
	return conv_polys
	
	
# ~ note: here, we take: k = n; 
def compute_trans_to_dHs_params(p, n, phi, ri_mat, neg_iri_mat, iM_dom, Zp, trans_v, ext_pol, R):  
	th1 = phi.powermod(n+1, p)		
	th2 = (th1 * phi) % p		
	omg = exact_conv_to_pmns__to_D1(th1, n, p, phi, ri_mat, neg_iri_mat, iM_dom, Zp, trans_v)			# a rep of 'phi^{n+1} mod p'
	dlt = exact_conv_to_pmns__to_D1(th2, n, p, phi, ri_mat, neg_iri_mat, iM_dom, Zp, trans_v)			# a rep of 'phi^{n+2} mod p'
	return (omg, dlt)

########################################################################

def compute_neg_inv_ri_mat(ri_mat, iM_dom):
	im = ri_mat.inverse()
	im = matrix(iM_dom, im)
	return -im

def compute_red_int_matrix(ri_poly, ext_pol, n, R):
	X = R.gen()
	ll = list(ri_poly)
	ll += [0]*(n-len(ll))
	bs = [ll]
	XiM = ri_poly
	for i in range(1,n):
		XiM = (X * XiM)%ext_pol
		ll = list(XiM)
		bs.append(ll + [0]*(n-len(ll)))
	return matrix(bs)

########################################################################
########################################################################

#important: assumes '0 <= v < base**n'
def compute_rep_in_base(v, base, n):
	v = Integer(v)
	rep = []
	while v != 0 :
		t = (v % base)
		rep.append(t)
		v //= base
	return rep + [0]*(n-len(rep))


#computes the lattice points in D_j
def compute_zero_reps_in_Dj(j, n, ri_mat):
	base = 2*j
	max_val = base**n
	trans_vect = vector(ZZ, [-j]*n)
	zero_points = []
	for v in range(max_val):
		vr = compute_rep_in_base(v, base, n)
		Z = vector(ZZ, vr) + trans_vect
		zero_points.append(Z*ri_mat)
	return zero_points


#computes the (2j)^n representations of 'v' in D_j 
#note: 'v' is taken modulo 'p' 
def compute_rep_set(v, p, n, omega_poly, ri_mat, neg_iri_mat, iM_dom, phi, ext_pol, trans_v, zero_trans_points, Zp, R):
	vr = exact_conv_to_pmns__to_D1(v, n, p, phi, ri_mat, neg_iri_mat, iM_dom, Zp, trans_v)
	vh = compute_rep_in_H(vr, n, omega_poly, ri_mat, neg_iri_mat, iM_dom, phi, ext_pol, trans_v, R)
	vr_set = [(vh+J) for J in zero_trans_points]
	return vr_set

########################################################################

#WARNING: should be the exact same function used for PMNS generation (see in the file 'pmns_generator.sage')
def build_lattice_base(p, n, gmm):
	b = []
	l = [p] + [0]*(n-1)
	b.append(l)
	m = identity_matrix(n)
	for i in range(1,n):
		t = (-gmm.powermod(i, p))%p
		if t%2 == 1 :
			t += p
		l = [t] + list(m[i][1:])
		b.append(l)
	bb = matrix(b)
	return bb.LLL(delta=1, algorithm='NTL:LLL')

########################################################################
########################################################################

def abs_matrix(mat):
	res = []
	for ll in mat:
		llabs = [abs(c) for c in ll]
		res.append(llabs)
	return matrix(res)


def compute_param_u(dom_j, p, n, delta, ext_pol_w, G):
	
	k = ceil(log(p,2)/n)
	bt = 1 << k
	c1 = n * (bt - 1)
	
	rho = dom_j*Integer(G.norm(1)) + 1
	
	f_add = (delta + 1)**2 
	c2 = ext_pol_w * f_add * (rho - 1)
	
	m = max(c1, c2)
	
	iG = G.inverse()
	u = ceil(m * (rho-1) * iG.norm(1))
	
	return u

########################################################################

def compute_extpol_w(n, ext_pol, R):
	x = R.gen()
	c = n-1
	tmp = x**n
	l1 = list(tmp%ext_pol)
	l1 = [abs(k) for k in l1]
	V = c * vector(l1 + [0]*(n-len(l1)))
	for d in range(n-2):
		c -= 1
		tmp *= x
		l1 = list(tmp%ext_pol)
		l1 = [abs(k) for k in l1]
		V += c * vector(l1 + [0]*(n-len(l1)))
	V += vector(range(1, n+1))
	return max(V)


def compute_translation_vect(n, u, ri_mat):
	iu = vector(ZZ, [-u]*n)
	return (iu * ri_mat)

########################################################################
########################################################################












