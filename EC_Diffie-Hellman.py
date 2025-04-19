###################################################################################################
#### IMPORTANT: we use Jacobian coordinates for short Weierstrass curves: y^2 = x^3 + a*x + b #####
###################################################################################################

#Note: pmns struct is [p, n, gmm, E, phi_log2, delta]

pmns = [76884956397045344220809746629001649093037950200943055203735601445031516197751, 5, 36181054844448274924615137064228865317857967132027946769785381863627930216959, [1, -1, 0, 0, 0, 1], 64, 5]  # a pmns for the brainpoolP256t1 modulus

#################### ELLIPTIC CURVE PARAMETERS #########################

#Equation: (E): y^2 = x^3 + ax + b, with:
Ea = int("7D5A0975FC2C3057EEF67530417AFFE7FB8055C126DC5C6CE94A4B44F330B5D9", 16)
Eb = int("26DC5C6CE94A4B44F330B5D9BBD77CBF958416295CF7E1CE6BCCDC18FF8C07B6", 16)

#Base point G(x,y) coordinates:
Gx = int("8BD2AEB9CB7E57CB2C4B482FFC81B7AFB9DE27E1E3BD23C23A4453BD9ACE3262", 16)
Gy = int("547EF835C3DAC4FD97F8461A14611DC9C27745132DED8E545C1D54C72F046997", 16)
std_affG = (Gx,Gy)

#G order: (it is also the curve order)
ord_G = int("A9FB57DBA1EEA9BC3E660A909D838D718C397AA3B561A6F7901E0E82974856A7", 16)

########################################################################

ZX.<x> = ZZ[] #or ZX.<X> = ZZ[]

load("ecdh-data.sage")	#WARNING: must always be loaded after 'pmns' and elliptic curve parameters initialisations as above

############### The Diffie-Hellman protocol execution ##################

PMNS_EC_Diffie_Hellman(base_G, ord_G, curve_A, curve_B, INFTY_P, ONE_rep, pmns_Mmat, pmns_iMmat, pmns_iM_dom, pmns_phi, pmns_E, pmns_Omega, pmns_phi2_rep, pmns_p, pmns_n, pmns_transV, ZX) 


########################################################################
##################### A series of checks ###############################
########################################################################

#making PMNS equality check to see if base_G is one the curve, with jacob (E): Y^2 = X^3 + aXZ^4 + bZ^6
print(PMNS_is_on_curve(base_G, curve_A, curve_B, pmns_Mmat, pmns_iMmat, pmns_iM_dom, pmns_phi, pmns_E, pmns_transV, ZX))

print(PMNS_is_on_curve(INFTY_P, curve_A, curve_B, pmns_Mmat, pmns_iMmat, pmns_iM_dom, pmns_phi, pmns_E, pmns_transV, ZX))

#checking backward conversion
std_backConv_G = from_PMNSJacPoint_to_stdAffPoint(base_G, pmns_n, pmns_gmm, pmns_phi, pmns_Zp, pmns_in_mont_domain, pmns_transV)
print(std_affG == std_backConv_G)

########################################################################
########################################################################

#checking curve arithmetic in the PMNS
GG = double_jacob_point(base_G, curve_A, pmns_Mmat, pmns_iMmat, pmns_iM_dom, pmns_phi, pmns_E, pmns_phi2_rep, pmns_transV, ZX)
GGG = add_jacob_points(base_G, GG, pmns_Mmat, pmns_iMmat, pmns_iM_dom, pmns_phi, pmns_E, pmns_phi2_rep, pmns_transV, ZX)
print(PMNS_is_on_curve(GG, curve_A, curve_B, pmns_Mmat, pmns_iMmat, pmns_iM_dom, pmns_phi, pmns_E, pmns_transV, ZX))
print(PMNS_is_on_curve(GGG, curve_A, curve_B, pmns_Mmat, pmns_iMmat, pmns_iM_dom, pmns_phi, pmns_E, pmns_transV, ZX))

#checking PMNS_ecsm_Montgomery_ladder:
sk = randint(1, ord_G-1)
kG = PMNS_ecsm_Montgomery_ladder(base_G, sk, curve_A, INFTY_P, pmns_Mmat, pmns_iMmat, pmns_iM_dom, pmns_phi, pmns_E, pmns_phi2_rep, pmns_transV, ZX)
print(PMNS_is_on_curve(kG, curve_A, curve_B, pmns_Mmat, pmns_iMmat, pmns_iM_dom, pmns_phi, pmns_E, pmns_transV, ZX))

########################################################################
########################################################################

#WARNING: wrong usages of point formulae: (give incorrect result) !!!!

add_jacob_points(base_G, INFTY_P, pmns_Mmat, pmns_iMmat, pmns_iM_dom, pmns_phi, pmns_E, pmns_phi2_rep, pmns_transV, ZX)
add_jacob_points(base_G, base_G, pmns_Mmat, pmns_iMmat, pmns_iM_dom, pmns_phi, pmns_E, pmns_phi2_rep, pmns_transV, ZX)
add_jacob_points(INFTY_P, INFTY_P, pmns_Mmat, pmns_iMmat, pmns_iM_dom, pmns_phi, pmns_E, pmns_phi2_rep, pmns_transV, ZX)
double_jacob_point(INFTY_P, curve_A, pmns_Mmat, pmns_iMmat, pmns_iM_dom, pmns_phi, pmns_E, pmns_phi2_rep, pmns_transV, ZX)

########################################################################

















