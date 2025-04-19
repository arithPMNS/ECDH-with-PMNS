load("diffie-hellman.sage")

########################################################################
############################# PMNS DATA ################################
########################################################################

pmns_p = pmns[0] 	 
pmns_n = pmns[1]
pmns_gmm = pmns[2]
pmns_E = ZX(pmns[3])
pmns_phi_log2 = pmns[4]  #must be >= 2 to ensure 'exact_conv_to_pmns__to_D1' consistency
pmns_dlt = pmns[5]

pmns_in_mont_domain = True	#must always be True here

pmns_Zp = GF(pmns_p)
pmns_convbase_log2 = ceil(log(pmns_p,2)/pmns_n)
pmns_phi = 1 << pmns_phi_log2
pmns_iM_dom = ZZ.quo(pmns_phi)
pmns_j = 1 
pmns_Mmat = build_lattice_base(pmns_p, pmns_n, pmns_gmm)  #WARNING: 'build_lattice_base' should be the exact same function used for PMNS generation.
pmns_Mmat_norm1 = Integer(pmns_Mmat.norm(1))
pmns_rho = 1 + pmns_Mmat_norm1
pmns_iMmat = compute_neg_inv_ri_mat(pmns_Mmat, pmns_iM_dom)

pmns_w = compute_extpol_w(pmns_n, pmns_E, ZX)
pmns_u = compute_param_u(pmns_j, pmns_p, pmns_n, pmns_dlt, pmns_w, pmns_Mmat)

pmns_transV = compute_translation_vect(pmns_n, pmns_u, pmns_Mmat)

pmns_conv_polys = compute_conv_polys__in_D1(pmns_p, pmns_n, pmns_convbase_log2, pmns_phi, pmns_Mmat, pmns_iMmat, pmns_iM_dom, pmns_Zp, pmns_transV, pmns_E, ZX, pmns_in_mont_domain)

(pmns_Omega, pmns_Delta) = compute_trans_to_dHs_params(pmns_p, pmns_n, pmns_phi, pmns_Mmat, pmns_iMmat, pmns_iM_dom, pmns_Zp, pmns_transV, pmns_E, ZX)

pmns_zero_reps = compute_zero_reps_in_Dj(pmns_j, pmns_n, pmns_Mmat)

pmns_phi2_rep = pmns_conv_polys[0] # a representation of 'pmns_phi^2' (since pmns_in_mont_domain = True); we will use it for exact coefficient reduction.

ONE_rep = conv_to_pmns__to_D1(1, pmns_n, pmns_convbase_log2, pmns_phi, pmns_conv_polys, pmns_Mmat, pmns_iMmat, pmns_iM_dom, pmns_Zp, pmns_transV)

ZERO_rep = pmns_zero_reps[randint(0,len(pmns_zero_reps)-1)]

# ~ print(f'Nb zeros: {len(pmns_zero_reps)}')
# ~ print(f'{pmns_zero_reps}')


########################################################################
############### ELLIPTIC CURVE PARAMETERS CONVERSION ###################
########################################################################

#IMORTANT: Conversion in THE MONTGOMERY DOMAIN is also done here

curve_A = conv_to_pmns__to_D1(Ea, pmns_n, pmns_convbase_log2, pmns_phi, pmns_conv_polys, pmns_Mmat, pmns_iMmat, pmns_iM_dom, pmns_Zp, pmns_transV)

curve_B = conv_to_pmns__to_D1(Eb, pmns_n, pmns_convbase_log2, pmns_phi, pmns_conv_polys, pmns_Mmat, pmns_iMmat, pmns_iM_dom, pmns_Zp, pmns_transV)

base_G = from_stdAffPoint_to_PMNSJacPoint(std_affG, pmns_n, pmns_convbase_log2, pmns_phi, pmns_conv_polys, pmns_Mmat, pmns_iMmat, pmns_iM_dom, pmns_Zp, pmns_transV)

#The point at infinity is the only point with a Z-coordinate equal to 0
#In Jacobian coordinates, it is represented by (1 : 1 : 0)
INFTY_P = (ONE_rep, ONE_rep, ZERO_rep) #remember conv in Mont domain is done

########################################################################
########################################################################

































