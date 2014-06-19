# coding=utf-8

# Adapted from Warren Moore's scripts.  The main function
# basic_info(curves) takes a list of elliptic curves and outputs a
# data file in the correct format.  Isogenous curves are computed and
# sorted.

# HNF comparison
def hnf_cmp(I, J):
        t = int(I.norm() - J.norm())
	if t:
		return t

	hnf_I = I.pari_hnf()
	hnf_J = J.pari_hnf()
	return int(hnf_I[1][0] - hnf_J[1][0])

# List of prime ideals
def prime_ideals(F, B):
	P = sum([p for p in [F.primes_above(p) for p in primes(B)]], [])
	P.sort(cmp = hnf_cmp)
	return P

# cached field data
Plists = {} # will be keyed by fields
Dlists = {} #
Glists = {} #
sigmas = {} #
labels = {} #
cm_j_invariants = {}
used_curves = {}

# Isogeny class comparison
ic_cmp = {}
def isog_class_cmp1(k, I, J):
	E_I = EllipticCurve(k, [int(ainvs[0]) + int(ainvs[1]) * k.gen() for ainvs in map(operator.methodcaller("split", ","), I[0].split()[6:11])])
	E_J = EllipticCurve(k, [int(ainvs[0]) + int(ainvs[1]) * k.gen() for ainvs in map(operator.methodcaller("split", ","), J[0].split()[6:11])])

	for p in Plists[k]:
		ap_I = ap(E_I, p)

		ap_J = ap(E_J, p)
		if ap_I != ap_J:
			return int(ap_I - ap_J)

	raise NotImplementedError("Bound on primes is too small to determine...")

def curve_cmp(E1,E2):
        ai1 = flatten([list(ai) for ai in E1.ainvs()])
        ai2 = flatten([list(ai) for ai in E2.ainvs()])
        return cmp(ai1,ai2)

def add_field(K):
    if K in used_curves:
        return
    Plists[K] = prime_ideals(K, 121)
    # Note: with bound only 100, curve [0,0,0,0,1] over Q(sqrt(-2))
    # passes the local isogeny test for l=499!
    D = K.discriminant()
    labels[K] = '2.0.%s.1' % str(-D) if D<0 else '2.2.%s.1' % str(D)
    Dlists[K] = D.abs()
    Glists[K] = K.galois_group()
    sigmas[K] = Glists[K][1]
# Get the CM j-invariants
    for d, f, j in cm_j_invariants_and_orders(K):
	cm_j_invariants[j] = d * (f ^ 2)
# Used curves
    used_curves[K] = {}
    ic_cmp[K] = lambda I,J: isog_class_cmp1(K,I,J)

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



# Check if we've already found this curve
def found(E, norm = None):
	if norm == None:
		norm = E.conductor().norm()
        K = E.base_field()
	if norm in used_curves[K]:
		for E2 in used_curves[K][norm]:
			if E.is_isomorphic(E2):
				return True
	return False

# Determine the isogeny class of E
isog_degrees = [2,3,5,7,11,13,17,19,23,29,31,41,47,59,71] # the "easy" l

# NB for curves with CM this will return a lot of primes which are not relevant!
def possible_isog_degrees(E, lmax=100):
        dlist = [0 if E.has_bad_reduction(P) else ap(E,P)^2-4*P.norm() for P in Plists[E.base_field()]]
        return [l for l in prime_range(lmax) if all([Mod(d,l).is_square() for d in dlist])]

import sys
def isog_class(E):
	isog_class = [E]
        # if E.j_invariant() in cm_j_invariants.keys():
        #         degs = isog_degrees
        # else:
        #         degs = possible_isog_degrees(E, lmax=1000)
        degs = possible_isog_degrees(E, lmax=1000)
        sys.stdout.write(" possible isogeny degrees: %s" % degs)
        bigdegs = [d for d in degs if d>100]
        if bigdegs:
                print " --Warning:  some large isogeny degrees will be tested!"
	isogenies = E.isogenies_prime_degree(degs)
        degs = list(set([phi.degree() for phi in isogenies]))
        print "... reduced to %s" % degs
	for isog in isogenies:
                E2 = isog.codomain()
                try:
                    E2 = E2.global_minimal_model()
                except:
                    pass

                if not any([E2.is_isomorphic(E3) for E3 in isog_class]):
			isog_class.append(E2)
			isogenies += E2.isogenies_prime_degree(degs)
        isog_class.sort(cmp=curve_cmp)
	return isog_class

def conj_curve(E,sigma):
        return EllipticCurve([sigma(a) for a in E.ainvs()])

def field_data(s):
    r"""
    Returns full field data from field label.
    """
    deg, r1, abs_disc, n = [int(c) for c in s.split(".")]
    sig = [r1, (deg-r1)//2]
    return [s, deg, sig, abs_disc]

def parse_NFelt(K, s):
    r"""
    Returns an element of K defined by the string s.
    """
    return K([QQ(c) for c in s.split(",")])

fields = {}
def field_from_label(lab):
        if lab in fields:
                return fields[lab]
        dummy, deg, sig, abs_disc = field_data(lab)
        d = ZZ(abs_disc)
        if sig[0]==0: d=-d
        x = polygen(QQ)
        t = d%4
        assert t in [0,1]
        pol = x^2 - t*x + (t-d)/4
        K = NumberField(pol, 'a')
        fields[lab] = K
        print "Created field from label %s: %s" % (lab,K)
        return K

def read_curves(infile):
        for L in file(infile).readlines():
                sys.stdout.write(L)
                data = L.split()
                if len(data)!=13:
                        print "line %s does not have 13 fields, skipping" % L
                        continue
                K = field_from_label(data[0])
                ainvs = [parse_NFelt(K,ai) for ai in data[6:11]]
                E = EllipticCurve(ainvs)
                yield E

# Basic info about the curves
def basic_info(curves, outfile = None):
        if outfile:
                outfile = file(outfile, mode="a")

	data = {}

	for E in curves:
                print "processing E = %s" % list(E.ainvs())
                k = E.base_field()
                add_field(k)
                D = Dlists[k]
                sigma = sigmas[k]
                used = used_curves[k]
                isog_class_cmp = ic_cmp[k]
                field_label = labels[k]
                if not k in data:
                        data[k] = {}
                data_k = data[k]

		# Get a global minimal model for E if possible
                try:
                        E = E.global_minimal_model()
                except:
                        pass
		N = E.conductor()
		norm = N.norm()

		if not found(E, norm):
			# Conductor
			hnf = N.pari_hnf()
			cond_label = "[%i,%s,%s]" % (norm, hnf[1][0], hnf[1][1])

			# Setup data
			if norm not in data_k:
				data_k[norm] = {}
			if hnf[1][0] not in data_k[norm]:
				data_k[norm][hnf[1][0]] = {}
			if hnf[1][1] not in data_k[norm][hnf[1][0]]:
				data_k[norm][hnf[1][0]][hnf[1][1]] = []
			else:
                                # This is only useful if we input a
                                # curve which is isogenous to one
                                # already processed but is not
                                # isomorphic to any previously seen,
                                # which only happens if the isog_class
                                # function produced an incomplete list
                                # from the earlier curve!
				for n, found_isog_class in enumerate(data_k[norm][hnf[1][0]][hnf[1][1]]):
					if E.is_isogenous(EllipticCurve(k, [int(ainvs[0]) + int(ainvs[1]) * k.gen() for ainvs in map(operator.methodcaller("split", ","), found_isog_class[0].split()[6:11])]), proof = False):
						ainvs = E.a_invariants()
						curve_data = found_isog_class[0].split()
						curve_data[3] = len(found_isog_class)+1
						curve_data[6:11] = ["%i,%i" % (ainvs[j][0], ainvs[j][1]) for j in xrange(0, 5)]
						data_k[norm][hnf[1][0]][hnf[1][1]][n].append(" ".join(curve_data))
						break

			# Let's find an isogeny class
			isogs = isog_class(E)
			if norm not in used:
				used[norm] = []
			used[norm] += isogs

                        # Q-curve? (isogeny class invariant)
                        if N != sigma(N):
                                q_curve = 0
                        elif all([ai[1]==0 for ai in ainvs]):
                                q_curve = 1
                        else:
                                Esigma = conj_curve(E,sigma)
                                q_curve = int(any([Esigma.is_isomorphic(E2) for E2 in isogs]))

			tmp = [] # list of output lines (with
                                 # placeholder for isog code, filled
                                 # in after sorting)

			for n, E2 in enumerate(isogs):
				# a-invs
				ainvs = E2.a_invariants()
				# Disc
				j = E2.j_invariant()
				disc = cm_j_invariants.get(j, 0)

				tmp.append("%s %s :isog %i %s %i %i,%i %i,%i %i,%i %i,%i %i,%i %i %s" % (field_label, cond_label, n + 1, cond_label, norm, ainvs[0][0], ainvs[0][1], ainvs[1][0], ainvs[1][1], ainvs[2][0], ainvs[2][1], ainvs[3][0], ainvs[3][1], ainvs[4][0], ainvs[4][1], disc, q_curve))
			data_k[norm][hnf[1][0]][hnf[1][1]].append(tmp)

	# Sort the isogeny classes
        ks = data.keys()
        print "fields: %s" % ks
        ks.sort()
        for k in ks:
            data_k = data[k]
            norms = data_k.keys()
            norms.sort()
            for norm in norms:
                data_k_n = data_k[norm]
		hnf0s = data_k_n.keys()
		hnf0s.sort()
		for hnf0 in hnf0s:
                        data_k_n_h = data_k_n[hnf0]
			hnf1s = data_k_n_h.keys()
			hnf1s.sort()
			for hnf1 in hnf1s:
                                dat = data_k_n_h[hnf1]
				dat.sort(cmp = isog_class_cmp)
				for n, isogs in enumerate(dat):
					isog_letter = chr(97 + n)
					for E_data in isogs:
                                                line = E_data.replace(":isog", isog_letter)
                                                if outfile:
                                                        outfile.write(line+'\n')
						print line
