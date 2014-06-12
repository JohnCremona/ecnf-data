# coding=utf-8

# HNF comparison
def hnf_cmp(I, J):
	if I.norm() != J.norm():
		return int(I.norm() - J.norm())
		
	hnf_I = I.pari_hnf()
	hnf_J = J.pari_hnf()
	
	return int(hnf_I[1][0] - hnf_J[1][0])
	
# List of prime ideals
def prime_ideals(F, B):
	P = sum([p for p in [F.primes_above(p) for p in primes(B)]], [])
	P.sort(cmp = hnf_cmp)
	return P

# Our field
#k.<i> = NumberField(x^2+1)
k=Qi
P = prime_ideals(k, 101)
D = abs(k.discriminant())
G = k.galois_group()
sigma = G[1]

# Used curves
used_curves = {}

# Get the CM j-invariants
cm_j_invariants = {}
for d, f, j in cm_j_invariants_and_orders(k):
	cm_j_invariants[j] = d * (f ^ 2)

# ap value
def ap(E, p):
	if E.has_good_reduction(p):
		return E.reduction(p).trace_of_frobenius()
	elif E.has_split_multiplicative_reduction(p):
		return 1
	elif E.has_nonsplit_multiplicative_reduction(p):
		return -1
	elif E.has_additive_reduction(p):
		return 0


# Isogeny class comparison
def isog_class_cmp(I, J):
	E_I = EllipticCurve(k, [int(ainvs[0]) + int(ainvs[1]) * k.gen() for ainvs in map(operator.methodcaller("split", ","), I[0].split()[6:11])])
	E_J = EllipticCurve(k, [int(ainvs[0]) + int(ainvs[1]) * k.gen() for ainvs in map(operator.methodcaller("split", ","), J[0].split()[6:11])])
	
	for p in P:
		ap_I = ap(E_I, p)
		ap_J = ap(E_J, p)
		if ap_I != ap_J:
			return int(ap_I - ap_J)
	
	raise NotImplementedError("Bound on primes is too small to determine...")

# Check if we've already found this curve
def found(E, norm = None):
	if norm == None:
		norm = E.conductor().norm()
	if norm in used_curves:
		for E2 in used_curves[norm]:
			if E.is_isomorphic(E2):
				return True
	return False

# Determine an isogeny class for E
def isog_class(E):
	isog_class = [E]
	isogenies = E.isogenies_prime_degree([2,3,5,7,11,13,17,19,23,29,31,41,47,59,71])
	for isog in isogenies:
		E2 = isog.codomain().global_minimal_model()
		flag = false
		for curve in isog_class:
			if E2.is_isomorphic(curve):
				flag = true
				break
		if not flag:
			isog_class.append(E2)
			isogenies += E2.isogenies_prime_degree([2,3,5,7,11,13,17,19,23,29,31,41,47,59,71])
	return isog_class

# Basic info about the curves
def basic_info(curves):
	data = {}
	
	for E in curves:
		# Get a global minimal model for E
		E = E.global_minimal_model()
		N = E.conductor()
		norm = N.norm()
		
		if not found(E, norm):
			# Conductor
			hnf = N.pari_hnf()
			cond_label = "[%i,%s,%s]" % (norm, hnf[1][0], hnf[1][1])
			
			# Setup data
			if norm not in data:
				data[norm] = {}
			if hnf[1][0] not in data[norm]:
				data[norm][hnf[1][0]] = {}
			if hnf[1][1] not in data[norm][hnf[1][0]]:
				data[norm][hnf[1][0]][hnf[1][1]] = []
			else:
				flag = False
				for n, found_isog_class in enumerate(data[norm][hnf[1][0]][hnf[1][1]]):
					if E.is_isogenous(EllipticCurve(k, [int(ainvs[0]) + int(ainvs[1]) * k.gen() for ainvs in map(operator.methodcaller("split", ","), found_isog_class[0].split()[6:11])]), proof = False):
						ainvs = E.a_invariants()
						curve_data = found_isog_class[0].split()
						curve_data[3] = len(found_isog_class)+1
						curve_data[6:11] = ["%i,%i" % (ainvs[j][0], ainvs[j][1]) for j in xrange(0, 5)]
						data[norm][hnf[1][0]][hnf[1][1]][n].append(" ".join(curve_data))
						flag = True
						break
				if flag == True:
					break
			
			# Let's find an isogeny class
			isogs = isog_class(E)
			if norm not in used_curves:
				used_curves[norm] = []
			used_curves[norm] += isogs
			
			tmp = []
			for n, E2 in enumerate(isogs):
				# a-invs
				ainvs = E2.a_invariants()
				
				# Disc
				j = E2.j_invariant()
				disc = cm_j_invariants[j] if j in cm_j_invariants else 0
				
				# Q-curve?
				q_curve = "?"
				if N != sigma(N):
					q_curve = 0
				elif ainvs[0][1] == 0 and ainvs[1][1] == 0 and ainvs[2][1] == 0 and ainvs[3][1] == 0 and ainvs[4][1] == 0:
					q_curve = 1
				else:
					try:
						q_curve = 1 if E2.is_isogenous(EllipticCurve(k, map(sigma, ainvs))) else 0
					except NotImplmentedError:
						pass
				
				tmp.append("2.0.%i.1 %s :isog %i %s %i %i,%i %i,%i %i,%i %i,%i %i,%i %i %s" % (D, cond_label, n + 1, cond_label, norm, ainvs[0][0], ainvs[0][1], ainvs[1][0], ainvs[1][1], ainvs[2][0], ainvs[2][1], ainvs[3][0], ainvs[3][1], ainvs[4][0], ainvs[4][1], disc, q_curve))
			data[norm][hnf[1][0]][hnf[1][1]].append(tmp)
	
	# Sort the isogeny classes
	norms = data.keys()
	norms.sort()
	for norm in norms:
		hnf0s = data[norm].keys()
		hnf0s.sort()
		for hnf0 in hnf0s:
			hnf1s = data[norm][hnf0].keys()
			hnf1s.sort()
			for hnf1 in hnf1s:
				data[norm][hnf0][hnf1].sort(cmp = isog_class_cmp)
				for n, isogs in enumerate(data[norm][hnf0][hnf1]):
					isog_letter = chr(97 + n)
					for E_data in isogs:
						print E_data.replace(":isog", isog_letter)
