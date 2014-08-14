# coding=utf-8
import sys
from sage.schemes.elliptic_curves.ell_curve_isogeny import fill_isogeny_matrix

# OLD version of curves.sage (up to August 2014).  Some of the
# functionality here for detecting CM and computing isogeny classes
# was implemented in Sage proper.  At time of writing (2014-08-14)
# this had not yet been accepted into Sage (latest version 6.3) but
# was awaiting review at http://trac.sagemath.org/ticket/16743 and the
# dependency http://trac.sagemath.org/ticket/16764.  The (shorter!)
# script curves.sage assumes the functionality provided by Sage-6.3 +
# those tickets.

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
labels = {} #
cm_j_invariants = {}
used_curves = {}

def curve_cmp(E1,E2):
        ai1 = flatten([list(ai) for ai in E1.ainvs()])
        ai2 = flatten([list(ai) for ai in E2.ainvs()])
        return cmp(ai1,ai2)

def add_field(K, charbound=200):
    if K in used_curves:
        return
    Plists[K] = prime_ideals(K, charbound)
    # Note: with bound only 100, curve [0,0,0,0,1] over Q(sqrt(-2))
    # passes the local isogeny test for l=499!
    # With bound 121 curve [0, 0, 0, -11, 14] over Q(sqrt(5)) tried l=617

    absD = K.discriminant().abs()
    s = K.signature()[0] # number of real places
    d = K.degree()
# Warning: number fields whose label's 4'th component is not 1 will
# not be handled correctly here
    labels[K] = '%s.%s.%s.1' % (str(d),str(s),str(absD))
    Dlists[K] = absD
    Glists[K] = K.galois_group(names='b')
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

# Functions which should be moved into the ell_number_field module:
def is_cm_j_invariant(j):
        if not j.is_integral():
                return False, None
        h = 1 if j in QQ else j.minpoly().degree()
        if h>100:
                raise NotImplementedError("CM data only available for class numbers up to 100")
        for d,f in cm_orders(h):
                pol = hilbert_class_polynomial(d*f^2)
                if pol(j)==0:
                        return True, (d,f)
        return False, None

def has_potential_cm(E):
        return is_cm_j_invariant(E.j_invariant())

def has_cm(E):
        flag, df = has_potential_cm(E)
        if not flag:
                return False, None
        d, f = df
        if E.base_field()(d).is_square():
                return True, d*f**2
        return False, None

# Determine the isogeny class of E
easy_isog_degrees = [2,3,5,7,11,13,17,19,23,29,31,41,47,59,71] # the "easy" l

def small_prime_value(Q):
        # returns a prime represented by a positive definite binary
        # quadratic form, not dividing its discriminant:
        d = Q.discriminant()
        for B in xsrange(10,1000,10):
                llist = list(Set([Q(x,y)
                                  for x in srange(-B,B) for y in srange(B)]))
                #llist = [l for l in llist if l.is_prime() and not l.divides(d)]
                llist = [l for l in llist if l.is_prime()]
                llist.sort()
                if llist:
                        return llist[0]
        raise ValueError("Unable to find a prime value of Q")

def possible_isog_degrees_CM(E, d):
        # First put in 2 and any odd primes l such that l* is a square:
        add_field(E.base_field())
        llist = K.absolute_discriminant().odd_part().prime_factors()
        llist = [l if l%4==1 else -l for l in llist]
        llist = [2] + [l.abs() for l in llist if K(l).is_square()]
        #print "first list: %s" % llist

        #Next put in primes which ramify in the CM order:
        llist += [l for l in d.prime_divisors() if not l in llist]
        #print "second list: %s" % llist
        llist = [2]
        # Now find primes (not dividing d) represented by each form of
        # discriminant d:
        Qs = BinaryQF_reduced_representatives(d, primitive_only=True)
        # discard principal form (q[0]=1) and one from an inverse pair:
        Qs = [q for q in Qs if q[0]>1 and q[1]>=0]
        for Q in Qs:
                l = small_prime_value(Q)
                if not l in llist:
                        llist.append(l)
        llist.sort()
        #print "final list: %s" % llist
        return llist


# NB for curves without CM we need to implement Billeray, but the
# following would only be incorrect if there is an isogeny of prime
# degree > 200 which is not very likely.

def possible_isog_degrees(E, lmax=200):
        add_field(E.base_field())
        K = E.base_field()
        flag, d = has_cm(E)
        if flag:
                return possible_isog_degrees_CM(E,d)
        return E.galois_representation().non_surjective()
        # dlist = [0 if E.has_bad_reduction(P) else ap(E,P)^2-4*P.norm()
        #          for P in Plists[E.base_field()]]
        # return [l for l in prime_range(lmax) if all([Mod(d,l).is_square()
        #                                                      for d in dlist])]

def isog_class(E, verbose=False):
        add_field(E.base_field())
        degs = possible_isog_degrees(E, lmax=1000)
        if verbose:
                sys.stdout.write(" possible isogeny degrees: %s" % degs)
                sys.stdout.flush()
        bigdegs = [d for d in degs if d>100]
        if bigdegs:
                print " --Warning:  some large isogeny degrees %s will be tested!" % degs
	isogenies = E.isogenies_prime_degree(degs)
        if verbose:
                sys.stdout.write(" -actual isogeny degrees: %s" % Set([phi.degree() for phi in isogenies]))
                sys.stdout.flush()
        # Add all new codomains to the list and collect degrees:
	curves = [E]
        ncurves = 1
        degs = []
        # triples (i,j,l) where curve i is l-isogenous to curve j
        triples = []
        def add_trip(t):
                for T in [t, [t[1],t[0],t[2]]]:
                        if not T in triples:
                                triples.append(T)
                                if verbose:
                                        sys.stdout.write(" -added triple %s..." % T)
                                        sys.stdout.flush()

        for phi in isogenies:
                E2 = phi.codomain()
                d = phi.degree()
                if not any([E2.is_isomorphic(E3) for E3 in curves]):
                        try:
                                E2 = E2.global_minimal_model()
                        except:
                                pass
                        curves.append(E2)
                        if verbose:
                                sys.stdout.write(" -added curve #%s (degree %s)..." % (ncurves,d))
                                sys.stdout.flush()
                        add_trip([0,ncurves,d])
                        ncurves += 1
                        if not d in degs:
                                degs.append(phi.degree())
        if verbose:
                sys.stdout.write("... relevant degrees: %s..." % degs)
                sys.stdout.write(" -now completing the isogeny class...")
                sys.stdout.flush()

        i = 1
        while i < ncurves:
                E1 = curves[i]
                if verbose:
                        sys.stdout.write(" -processing curve #%s..." % i)
                        sys.stdout.flush()

                isogenies = E1.isogenies_prime_degree(degs)

                for phi in isogenies:
                        E2 = phi.codomain()
                        d = phi.degree()
                        js = [j for j,E3 in enumerate(curves) if E2.is_isomorphic(E3)]
                        if js: # seen codomain already
                                add_trip([i,js[0],d])
                        else:
                                try:
                                        E2 = E2.global_minimal_model()
                                except:
                                        pass
                                curves.append(E2)
                                if verbose:
                                        sys.stdout.write(" -added curve #%s..." % ncurves)
                                        sys.stdout.flush()
                                add_trip([i,ncurves,d])
                                ncurves += 1
                i += 1

        scurves = sorted(curves,cmp=curve_cmp)
        perm = dict([(i,scurves.index(E)) for i,E in enumerate(curves)])
        if verbose:
                print "Sorting permutation = %s" % perm

        mat = MatrixSpace(ZZ,ncurves)(0)
        for i,j,l in triples:
                mat[perm[i],perm[j]] = l
        mat = fill_isogeny_matrix(mat)

        if verbose:
                print
                print("... isogeny class has size %s" % ncurves)
                #print("Matrix = \n%s" % mat)
	return curves, mat

def is_Galois_invariant(N):
        try:
                K = N.number_field()
        except AttributeError:
                try:
                        K = N.parent()
                except AttributeError:
                        raise ValueError("unable to determine field from %s" % N)
        if K is QQ: return True
        add_field(K)
        G = Glists[K]
        NL = G[0](N) # base-change to Galois closure
        return all([sigma(N)==NL for sigma in G.gens()])

def conj_curve(E,sigma):
        return EllipticCurve([sigma(a) for a in E.ainvs()])

def is_Q_curve(E,isogs=None):
        K = E.base_field()
        if K is QQ: return True
        if not is_Galois_invariant(E.conductor()):
                return False
        if all([is_Galois_invariant(a) for a in E.ainvs()]):
                return True
        add_field(K)
        G = Glists[K]
        EL = conj_curve(E,G[0]) # base-change to Galois closure
        if isogs is None:
                isogs, mat = isog_class(E)
        return all([any([conj_curve(E,sigma).is_isomorphic(EL) for E2 in isogs]) for sigma in G.gens()])

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

def ainvs_to_string(ainvs):
        return " ".join([",".join([str(c) for c in list(ai)]) for ai in ainvs])

def ainvs_from_strings(K, ainv_string_list):
        return [parse_NFelt(K,ai) for ai in ainv_string_list]

def curve_from_strings(K, ainv_string_list):
        return EllipticCurve(ainvs_from_strings(K,ainv_string_list))

# Isogeny class comparison
ic_cmp = {}
def isog_class_cmp1(k, I, J):
	E_I = curve_from_strings(k,I[0].split()[6:11])
	E_J = curve_from_strings(k,J[0].split()[6:11])

	for p in Plists[k]:
		c = int(ap(E_I, p) - ap(E_J, p))
		if c: return c

	raise NotImplementedError("Bound on primes is too small to determine...")


fields = {}
def field_from_label(lab):
        if lab in fields:
                return fields[lab]
        dummy, deg, sig, abs_disc = field_data(lab)
        x = polygen(QQ)
        if deg==2:
                d = ZZ(abs_disc)
                if sig[0]==0: d=-d
                t = d%4
                assert t in [0,1]
                pol = x^2 - t*x + (t-d)/4
        elif lab=='3.1.23.1':
                pol = x**3 - x**2 +1
        else:
                raise NotImplementedError("cannot yet handle field %s" % lab)
        K = NumberField(pol, 'a')
        fields[lab] = K
        print "Created field from label %s: %s" % (lab,K)
        return K

def read_curves(infile):
        for L in file(infile).readlines():
                #sys.stdout.write(L)
                data = L.split()
                if len(data)!=13:
                        print "line %s does not have 13 fields, skipping" % L
                        continue
                K = field_from_label(data[0])
                E = curve_from_strings(K, data[6:11])
                yield E

# Basic info about the curves
def basic_info(curves, outfile = None, verbose=0):
        if outfile:
                outfile = file(outfile, mode="a")

	data = {}

	for E in curves:
                if verbose>0:
                        print("processing E = %s..." % list(E.ainvs()))
                k = E.base_field()
                add_field(k)
                D = Dlists[k]
                G = Glists[k]
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

		if found(E, norm):
                        if verbose>0:
                                print(" -- isogenous to a previous curve")
                else:
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
                                ainvs = E.a_invariants()
				for n, found_isog_class in enumerate(data_k[norm][hnf[1][0]][hnf[1][1]]):
                                        curve_data = found_isog_class[0].split()
					if E.is_isogenous(curve_from_strings(k, curve_data[6:11]), proof = False):
						curve_data[3] = len(found_isog_class)+1
						curve_data[6:11] = [",".join([str(c) for c in ai]) for ai in ainvs]
						data_k[norm][hnf[1][0]][hnf[1][1]][n].append(" ".join(curve_data))
						break

			# Let's find an isogeny class
			isogs, mat = isog_class(E, verbose>1)
			if norm not in used:
				used[norm] = []
			used[norm] += isogs

                        # Q-curve? (isogeny class invariant)
                        q_curve = int(is_Q_curve(E))

			tmp = [] # list of output lines (with
                                 # placeholder for isog code, filled
                                 # in after sorting)

			for n, E2 in enumerate(isogs):
				# a-invs
				ainvs = E2.a_invariants()
                                ainv_string = ainvs_to_string(ainvs)
				# Disc
				j = E2.j_invariant()
				disc = cm_j_invariants.get(j, 0)

				tmp.append("%s %s :isog %i %s %i %s %i %i" % (field_label, cond_label, n + 1, cond_label, norm, ainv_string, disc, q_curve))
                        #print "appending %s curves" % len(tmp)
			data_k[norm][hnf[1][0]][hnf[1][1]].append(tmp)

	# Sort the isogeny classes
        ks = data.keys()
        if verbose>0:
                print
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
						if verbose>0:
                                                        print line
