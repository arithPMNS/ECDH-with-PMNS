load("pmns_arith_ops.sage")

###################################################################################################
#### IMPORTANT: we use Jacobian coordinates for short Weierstrass curves: y^2 = x^3 + a*x + b #####
###################################################################################################

# ~ assumes that P has at least two coordinates
def print_ec_point(P):
	rr = '(' + str(P[0]) + ',\n'
	for c in P[1:-1]:
		rr += ' ' + str(c) + ',\n'
	rr += ' ' + str(P[-1]) + ')'
	print(rr)


def print_ec_point_polys(P, ZX):
	for el in P:
		print(f'  {ZX(list(el))}')
	

# ~ assumes Jacobian coordinates for the base point G, with coord. already in the PMNS
# ~ the share secret is the X affine coordinate in the fundamental region H
def PMNS_EC_Diffie_Hellman(G, ord_G, curve_A, curve_B, infty_P, ONE_rep, ri_mat, neg_iri_mat, iM_dom, phi, ext_pol, omega_poly, phi2_rep, p, n, trans_v, ZX):

	#-------------------------------------------------------------------
	print('\n-------------------------------------------------------------------')
	
	print(f'Base point G is on curve: {PMNS_is_on_curve(G, curve_A, curve_B, ri_mat, neg_iri_mat, iM_dom, phi, ext_pol, trans_v, ZX)}')
		
	#-------------------------------------------------------------------
	print('-------------------------------------------------------------------')
		
	alice_sk = randint(1, ord_G-1) 
	alice_kG = PMNS_ecsm_Montgomery_ladder(G, alice_sk, curve_A, infty_P, ri_mat, neg_iri_mat, iM_dom, phi, ext_pol, phi2_rep, trans_v, ZX)
	print(f'\nAlice secret scalar: \n{alice_sk}')
	print('Alice public point coordinates:')
	print_ec_point_polys(alice_kG, ZX)
	# ~ print_ec_point(alice_kG)
	#Note: Alice sends 'alice_kG' to Bob
	
	bob_sk = randint(1, ord_G-1) 
	bob_kG = PMNS_ecsm_Montgomery_ladder(G, bob_sk, curve_A, infty_P, ri_mat, neg_iri_mat, iM_dom, phi, ext_pol, phi2_rep, trans_v, ZX)
	print(f'\nBob secret scalar: \n{bob_sk}')
	print('Bob public point coordinates:')
	print_ec_point_polys(bob_kG, ZX)
	# ~ print_ec_point(bob_kG)
	#Note: Bob sends 'bob_kG' to Alice
	
	print('\nChecking if public points are on the curve:')
	print(f'  Alice public point: {PMNS_is_on_curve(alice_kG, curve_A, curve_B, ri_mat, neg_iri_mat, iM_dom, phi, ext_pol, trans_v, ZX)}')
	print(f'  Bob public point  : {PMNS_is_on_curve(bob_kG, curve_A, curve_B, ri_mat, neg_iri_mat, iM_dom, phi, ext_pol, trans_v, ZX)}')
	
	#-------------------------------------------------------------------
	print('\n-------------------------------------------------------------------\n')
	
	#Alice computes:
	alice_secret = PMNS_ecsm_Montgomery_ladder(bob_kG, alice_sk, curve_A, infty_P, ri_mat, neg_iri_mat, iM_dom, phi, ext_pol, phi2_rep, trans_v, ZX)
	alice_aff = from_PMNSJacPoint_to_PMNSAffPoint(alice_secret, p, ONE_rep, ri_mat, neg_iri_mat, iM_dom, phi, ext_pol, trans_v, ZX)
	print('Alice secret point coordinates:')
	print_ec_point_polys(alice_aff, ZX)
	# ~ print_ec_point(alice_aff)
	
	#Bob computes:
	bob_secret = PMNS_ecsm_Montgomery_ladder(alice_kG, bob_sk, curve_A, infty_P, ri_mat, neg_iri_mat, iM_dom, phi, ext_pol, phi2_rep, trans_v, ZX)
	bob_aff = from_PMNSJacPoint_to_PMNSAffPoint(bob_secret, p, ONE_rep, ri_mat, neg_iri_mat, iM_dom, phi, ext_pol, trans_v, ZX)
	print('\nBob secret point coordinates:')
	print_ec_point_polys(bob_aff, ZX)
	# ~ print_ec_point(bob_aff)
	
	#-------------------------------------------------------------------
	print('\n-------------------------------------------------------------------')
	
	#Conversion of the affine coordinates to the fundamental region H, where the representations are unique
	
	alice_affX_in_H = compute_rep_in_H(alice_aff[0], n, omega_poly, ri_mat, neg_iri_mat, iM_dom, phi, ext_pol, trans_v, ZX)
	print(f'\nAlice secret affX in H: \n{alice_affX_in_H}')
	
	bob_affX_in_H = compute_rep_in_H(bob_aff[0], n, omega_poly, ri_mat, neg_iri_mat, iM_dom, phi, ext_pol, trans_v, ZX)
	print(f'\nBob secret affX in H: \n{bob_affX_in_H}')
	
	print(f'\nAre same polys in H? {alice_affX_in_H == bob_affX_in_H}')

	#-------------------------------------------------------------------
	print('\n-------------------------------------------------------------------\n')
	
	return;

########################################################################

# ~ assumes Jacobian coordinates for P, with coord. already in the PMNS (and in Mont domain)
# ~ also assumes curve_A and curve_B in the PMNS (and in Mont domain)
# ~ Note: Jacobian coordinates equation is Y^2 = X^3 + aXZ^4 + bZ^6
def PMNS_is_on_curve(P, curve_A, curve_B, ri_mat, neg_iri_mat, iM_dom, phi, ext_pol, trans_v, ZX):
	
	X1 = P[0]
	Y1 = P[1]
	Z1 = P[2]
	
	Y2 = pmns_trans_mont_mult(Y1, Y1, ri_mat, neg_iri_mat, iM_dom, phi, ext_pol, trans_v, ZX)
	
	X2 = pmns_trans_mont_mult(X1, X1, ri_mat, neg_iri_mat, iM_dom, phi, ext_pol, trans_v, ZX)
	X3 = pmns_trans_mont_mult(X2, X1, ri_mat, neg_iri_mat, iM_dom, phi, ext_pol, trans_v, ZX)
	
	Z2 = pmns_trans_mont_mult(Z1, Z1, ri_mat, neg_iri_mat, iM_dom, phi, ext_pol, trans_v, ZX)
	Z4 = pmns_trans_mont_mult(Z2, Z2, ri_mat, neg_iri_mat, iM_dom, phi, ext_pol, trans_v, ZX)
	Z6 = pmns_trans_mont_mult(Z4, Z2, ri_mat, neg_iri_mat, iM_dom, phi, ext_pol, trans_v, ZX)
	
	tmpA = pmns_trans_mont_mult(curve_A, X1, ri_mat, neg_iri_mat, iM_dom, phi, ext_pol, trans_v, ZX)
	tmpA = pmns_trans_mont_mult(tmpA, Z4, ri_mat, neg_iri_mat, iM_dom, phi, ext_pol, trans_v, ZX)
	
	tmpB = pmns_trans_mont_mult(curve_B, Z6, ri_mat, neg_iri_mat, iM_dom, phi, ext_pol, trans_v, ZX)
	
	tmpAB = X3 + tmpA + tmpB
	
	return equality_check(Y2, tmpAB, ri_mat, neg_iri_mat, iM_dom, phi, trans_v)

########################################################################

# ~ assumes affine coordinates (x,y) for P
# ~ note: conv to Montgomrey domain is also done for the coordinates
def from_stdAffPoint_to_PMNSJacPoint(P, n, convbase_log2, phi, conv_polys, ri_mat, neg_iri_mat, iM_dom, Zp, trans_v):
	
	one_rep = conv_to_pmns__to_D1(1, n, convbase_log2, phi, conv_polys, ri_mat, neg_iri_mat, iM_dom, Zp, trans_v)
	
	jPx = conv_to_pmns__to_D1(P[0], n, convbase_log2, phi, conv_polys, ri_mat, neg_iri_mat, iM_dom, Zp, trans_v)
	jPy = conv_to_pmns__to_D1(P[1], n, convbase_log2, phi, conv_polys, ri_mat, neg_iri_mat, iM_dom, Zp, trans_v)
	
	jP = (jPx, jPy, one_rep)
	
	return jP


def from_PMNSJacPoint_to_stdAffPoint(P, n, gmm, phi, Zq, in_mont_domain):
	
	std_X = conv_from_pmns(P[0], n, gmm, phi, Zq, in_mont_domain)
	std_Y = conv_from_pmns(P[1], n, gmm, phi, Zq, in_mont_domain)
	std_Z = conv_from_pmns(P[2], n, gmm, phi, Zq, in_mont_domain)
	Z_inv = Zq(1/std_Z)
	Z_inv2 = Zq(Z_inv**2)
	Z_inv3 = Zq(Z_inv2 * Z_inv)
	aff_X = Zq(std_X * Z_inv2)
	aff_Y = Zq(std_Y * Z_inv3)
	return (aff_X, aff_Y)
	

def from_PMNSJacPoint_to_PMNSAffPoint(P, p, ONE_rep, ri_mat, neg_iri_mat, iM_dom, phi, ext_pol, trans_v, ZX):
	
	Z = P[2]
	Z_inv = PMNS_modInverse(Z, p, ONE_rep, ri_mat, neg_iri_mat, iM_dom, phi, ext_pol, trans_v, ZX)
	Zinv2 = pmns_trans_mont_mult(Z_inv, Z_inv, ri_mat, neg_iri_mat, iM_dom, phi, ext_pol, trans_v, ZX)
	Zinv3 = pmns_trans_mont_mult(Zinv2, Z_inv, ri_mat, neg_iri_mat, iM_dom, phi, ext_pol, trans_v, ZX)
	aff_X = pmns_trans_mont_mult(P[0], Zinv2, ri_mat, neg_iri_mat, iM_dom, phi, ext_pol, trans_v, ZX)
	aff_Y = pmns_trans_mont_mult(P[1], Zinv3, ri_mat, neg_iri_mat, iM_dom, phi, ext_pol, trans_v, ZX)
	return (aff_X, aff_Y)
	
	
def from_PMNSAffPoint_to_stdAffPoint(P, n, gmm, phi, Zq, in_mont_domain):
	
	std_x = conv_from_pmns(P[0], n, gmm, phi, Zq, in_mont_domain)
	std_y = conv_from_pmns(P[1], n, gmm, phi, Zq, in_mont_domain)
	return (std_x, std_y)

########################################################################

# ~ assumes Jacobian coordinates for P and inf_P, with coord. already in the PMNS
# ~ IMPORTANT: make sure P is NOT the point at infinity
def PMNS_ecsm_Montgomery_ladder(P, k, curve_A, infty_P, ri_mat, neg_iri_mat, iM_dom, phi, ext_pol, phi2_rep, trans_v, ZX):
	
	k_bits = [int(b) for b in bin(k)[2:]] #bits from left to right
	#IMPORTANT: For correct use of point formulae, the MSB of k is assumed to be 1 and is integrated in RR initialisation:
	PP = double_jacob_point(P, curve_A, ri_mat, neg_iri_mat, iM_dom, phi, ext_pol, phi2_rep, trans_v, ZX)
	RR = [P, PP]
	for b in k_bits[1:]: #first bit (which is 1) already used above
		RR[1-b] = add_jacob_points(RR[1-b], RR[b], ri_mat, neg_iri_mat, iM_dom, phi, ext_pol, phi2_rep, trans_v, ZX)
		RR[b] = double_jacob_point(RR[b], curve_A, ri_mat, neg_iri_mat, iM_dom, phi, ext_pol, phi2_rep, trans_v, ZX)
	return RR[0]
	

# ~ assumes A already in the PMNS
# ~ also assumes A and ONE_rep in the Montgomrey domain; so ONE_rep represents 'phi' 
def PMNS_Montgomery_ladder(A, k, ONE_rep, ri_mat, neg_iri_mat, iM_dom, phi, ext_pol, trans_v, ZX):
	
	k_bits = [int(b) for b in bin(k)[2:]] #bits from left to right
	RR = [ONE_rep, A]
	for b in k_bits:
		RR[1-b] = pmns_trans_mont_mult(RR[1-b], RR[b], ri_mat, neg_iri_mat, iM_dom, phi, ext_pol, trans_v, ZX)
		RR[b] = pmns_trans_mont_mult(RR[b], RR[b], ri_mat, neg_iri_mat, iM_dom, phi, ext_pol, trans_v, ZX)
	return RR[0]
	

# ~ assumes that A is in the Montgomrey domain and is NOT a representation of 0
# ~ Note: we use Euler's theorem
def PMNS_modInverse(A, p, ONE_rep, ri_mat, neg_iri_mat, iM_dom, phi, ext_pol, trans_v, ZX):
	k = p - 2 
	A_inv = PMNS_Montgomery_ladder(A, k, ONE_rep, ri_mat, neg_iri_mat, iM_dom, phi, ext_pol, trans_v, ZX)
	return A_inv

########################################################################

# ~ assumes Jacobian coordinates for P and Q, with coord. already in the PMNS
def add_jacob_points(P, Q, ri_mat, neg_iri_mat, iM_dom, phi, ext_pol, phi2_rep, trans_v, ZX): 
	
	X1 = P[0]
	Y1 = P[1]
	Z1 = P[2]
	X2 = Q[0]
	Y2 = Q[1]
	Z2 = Q[2]
	
	Z1Z1 = pmns_trans_mont_mult(Z1, Z1, ri_mat, neg_iri_mat, iM_dom, phi, ext_pol, trans_v, ZX)
	Z2Z2 = pmns_trans_mont_mult(Z2, Z2, ri_mat, neg_iri_mat, iM_dom, phi, ext_pol, trans_v, ZX)
	
	U1 = pmns_trans_mont_mult(X1, Z2Z2, ri_mat, neg_iri_mat, iM_dom, phi, ext_pol, trans_v, ZX)
	U2 = pmns_trans_mont_mult(X2, Z1Z1, ri_mat, neg_iri_mat, iM_dom, phi, ext_pol, trans_v, ZX)
	
	S1 = pmns_trans_mont_mult(Z2, Z2Z2, ri_mat, neg_iri_mat, iM_dom, phi, ext_pol, trans_v, ZX)
	S1 = pmns_trans_mont_mult(S1, Y1, ri_mat, neg_iri_mat, iM_dom, phi, ext_pol, trans_v, ZX)
	
	S2 = pmns_trans_mont_mult(Z1, Z1Z1, ri_mat, neg_iri_mat, iM_dom, phi, ext_pol, trans_v, ZX)
	S2 = pmns_trans_mont_mult(S2, Y2, ri_mat, neg_iri_mat, iM_dom, phi, ext_pol, trans_v, ZX)

	H = U2 - U1
	I = 2 * H
	I = pmns_trans_mont_mult(I, I, ri_mat, neg_iri_mat, iM_dom, phi, ext_pol, trans_v, ZX)
	J = pmns_trans_mont_mult(H, I, ri_mat, neg_iri_mat, iM_dom, phi, ext_pol, trans_v, ZX)
	V = pmns_trans_mont_mult(U1, I, ri_mat, neg_iri_mat, iM_dom, phi, ext_pol, trans_v, ZX)
	
	r = 2 * (S2 - S1)
	r2 = pmns_trans_mont_mult(r, r, ri_mat, neg_iri_mat, iM_dom, phi, ext_pol, trans_v, ZX)
	
	X3 = r2 - J - 2*V
	X3 = exact_red_coeffs(X3, ri_mat, neg_iri_mat, iM_dom, phi, ext_pol, phi2_rep, trans_v, ZX)
	
	tmp1 = pmns_trans_mont_mult(r, V-X3, ri_mat, neg_iri_mat, iM_dom, phi, ext_pol, trans_v, ZX)
	tmp2 = pmns_trans_mont_mult(S1, J, ri_mat, neg_iri_mat, iM_dom, phi, ext_pol, trans_v, ZX)
	Y3 = tmp1 - 2*tmp2
	Y3 = exact_red_coeffs(Y3, ri_mat, neg_iri_mat, iM_dom, phi, ext_pol, phi2_rep, trans_v, ZX)
	
	tmp3 = pmns_trans_mont_mult(Z1+Z2, Z1+Z2, ri_mat, neg_iri_mat, iM_dom, phi, ext_pol, trans_v, ZX)
	tmp4 = tmp3 - Z1Z1 - Z2Z2
	Z3 = pmns_trans_mont_mult(tmp4, H, ri_mat, neg_iri_mat, iM_dom, phi, ext_pol, trans_v, ZX)
	
	RR = (X3, Y3, Z3)
	
	return RR


# ~ assumes Jacobian coordinates for P, with coord. already in the PMNS
# ~ also assumes curve_A in the Montgomrey domain in the PMNS
def double_jacob_point(P, curve_A, ri_mat, neg_iri_mat, iM_dom, phi, ext_pol, phi2_rep, trans_v, ZX): 
	
	X1 = P[0]
	Y1 = P[1]
	Z1 = P[2]
	 
	XX = pmns_trans_mont_mult(X1, X1, ri_mat, neg_iri_mat, iM_dom, phi, ext_pol, trans_v, ZX)
	YY = pmns_trans_mont_mult(Y1, Y1, ri_mat, neg_iri_mat, iM_dom, phi, ext_pol, trans_v, ZX)
	YYYY = pmns_trans_mont_mult(YY, YY, ri_mat, neg_iri_mat, iM_dom, phi, ext_pol, trans_v, ZX)
	ZZ = pmns_trans_mont_mult(Z1, Z1, ri_mat, neg_iri_mat, iM_dom, phi, ext_pol, trans_v, ZX)
	
	tmp1 = pmns_trans_mont_mult(X1 + YY, X1 + YY, ri_mat, neg_iri_mat, iM_dom, phi, ext_pol, trans_v, ZX)
	S = 2 * (tmp1 - XX - YYYY)
	
	ZZZZ = pmns_trans_mont_mult(ZZ, ZZ, ri_mat, neg_iri_mat, iM_dom, phi, ext_pol, trans_v, ZX)
	tmp_a = pmns_trans_mont_mult(curve_A, ZZZZ, ri_mat, neg_iri_mat, iM_dom, phi, ext_pol, trans_v, ZX)
	M = 3*XX + tmp_a
	
	MM = pmns_trans_mont_mult(M, M, ri_mat, neg_iri_mat, iM_dom, phi, ext_pol, trans_v, ZX)
	T = MM - 2*S
	
	X3 = exact_red_coeffs(T, ri_mat, neg_iri_mat, iM_dom, phi, ext_pol, phi2_rep, trans_v, ZX)
	
	tmp2 = pmns_trans_mont_mult(M, S-T, ri_mat, neg_iri_mat, iM_dom, phi, ext_pol, trans_v, ZX)
	Y3 = tmp2 - 8*YYYY
	Y3 = exact_red_coeffs(Y3, ri_mat, neg_iri_mat, iM_dom, phi, ext_pol, phi2_rep, trans_v, ZX)
	
	tmp3 = pmns_trans_mont_mult(Y1+Z1, Y1+Z1, ri_mat, neg_iri_mat, iM_dom, phi, ext_pol, trans_v, ZX)
	Z3 = tmp3 - YY - ZZ
	Z3 = exact_red_coeffs(Z3, ri_mat, neg_iri_mat, iM_dom, phi, ext_pol, phi2_rep, trans_v, ZX)
	
	RR = (X3, Y3, Z3)
	
	return RR

########################################################################






