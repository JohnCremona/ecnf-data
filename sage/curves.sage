# coding=utf-8
import sys
from sage.databases.cremona import cremona_letter_code

# Adapted from Warren Moore's scripts.  The main function
# basic_info(curves) takes a list of elliptic curves and outputs a
# data file in the correct format.  Isogenous curves are computed and
# sorted.

# NB requires functionality of Sage-6.3 and code for computing isogeny
# classes as implemented in Sage, but at present (2014-08-14) this has
# not yet all been accepted into Sage (latest version 6.3) but is
# awaiting review at http://trac.sagemath.org/ticket/16743 and the
# dependency http://trac.sagemath.org/ticket/16764.

# cached field data: these will all be keyed by fields k, values as shown:
Plists = {} # list of primes of k, sorted by norm
Dlists = {} # |disc(k)|
Glists = {} # Gal(k/Q)
labels = {} # label of k
nf_data = {} # newform data for k (if available)
used_curves = {} # dict with key by norm(conductor), value a list of
                 # curves so far processed witth that norm-conductor
ic_cmp = {} # dict whose values is the isogeny class cmp function
cm_counts = {} # dict with keys conductor labels, values counts of
               # classes with rational CM (imaginary quadratic fields
               # only) for labelling of these, as they do not have
               # associated Bianchi newforms
#
cm_j_invariants = {} # key: CM j-invariants (in any field)
                     # value: associated negative discriminant

def add_field(K, field_label=None, prime_norm_bound=200):
    if K in used_curves:
        return

    PP = K.primes_of_bounded_norm(prime_norm_bound)
    PP.sort(key=lambda P: (P.norm(),P.pari_hnf().sage()))
    Plists[K] = PP

    absD = K.discriminant().abs()
    s = K.signature()[0] # number of real places
    d = K.degree()
# Warning: number fields whose label's 4'th component is not 1 will
# not be handled correctly here; this is only an issue when there is
# more than one field of given signature and abs(disc), so fine for
# quadratic fields.  This problem first hit for field 4.4.16448.2.
# When we processed te curves for that field they were in a different
# input file so were not mixed up with curves for 4.4.16448.1 luckily,
# and manual editing of the output was sufficient.

    if field_label==None:
        field_label = '%s.%s.%s.1' % (str(d),str(s),str(absD))
    labels[K] = field_label
    Dlists[K] = absD
    Glists[K] = K.galois_group(names='b')
    for dd, f, j in cm_j_invariants_and_orders(K):
	cm_j_invariants[j] = dd * (f ^ 2)
    used_curves[K] = {}
    ic_cmp[K] = lambda I,J: isog_class_cmp1(K,I,J)
    nf_data[K] = None
    cm_counts[K] = {}
    if d==2 and s==0 and absD in [3,4,7,8,11]:
            from nfscripts import read_newform_data, nf_filename
            nf_data[K] = read_newform_data(nf_filename(absD))

def ap(E, p):
        r"""
        Return a_p(E).

        INPUT:

        - ``E`` - an elliptic curve defined over a number field `k`;

        - ``p`` - a prime odeal of `k`.

        OUTPUT:

        `a_p(E)`: the trace of Frobenius of `E` at `p` if `E` has good
        reduction, otherwise the appropriate L-series coefficient
        depending on the type of bad reduction.
        """
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
        r"""
        Return ``True`` iff this isomorphism class has been encountered already.
        """
	if norm is None:
		norm = E.conductor().norm()
        K = E.base_field()
	if norm in used_curves[K]:
		for E2 in used_curves[K][norm]:
			if E.is_isomorphic(E2):
				return E2
	return None

# Functions for testing if E is a Q-curve

def is_Galois_invariant(N, field_label=None):
        r"""
        Return ``True`` if this number field element or ideal is Galois-invariant.
        """
        try:
                K = N.number_field()
        except AttributeError:
                try:
                        K = N.parent()
                except AttributeError:
                        raise ValueError("unable to determine field from %s" % N)
        if K is QQ: return True
        add_field(K, field_label=field_label)
        G = Glists[K]
        NL = G[0](N) # base-change to Galois closure
        return all([sigma(N)==NL for sigma in G.gens()])

def conj_curve(E,sigma):
        r"""
        Return the Galois conjugate elliptic curve under sigma.
        """
        return EllipticCurve([sigma(a) for a in E.ainvs()])

def is_Q_curve(E, field_label=None):
        r"""
        Return True if this elliptic curve is isogenous to all its
        Galois conjugates.

        Note: if the base field K is not Galois we compute the Galois
        group of its Galois closure L, and test for isogeny over L.
        Is this right?  Following Elkies ('On elliptic K-curves',
        2004) this is the correct set of Galois conjugates but we
        should be looking for isogeny over the algebraic closure.

        If E does not have CM (defined over L) and there is an isogeny
        phi from E to E' where E' is defined over L but phi need not
        be, then considering the composite of phi with the dual of its
        Galois conjugates shows that each of these conjugates is equal
        to phi up to sign.  If the signs are all +1 then phi is also
        defined over L, but otherwise it is only defined over a
        quadratic extension M of L.  In that case, replacing E' by its
        quadratic twist and phi by its composite with the isomorphism
        (defined over M) from E' to its twist gives a new curve
        defined over L and isogenous to E via an L-rational isogeny of
        the same degree as the original phi.  For our test for being a
        Q-curve to be correct (in this non-CM case) -- i.e., agree
        with Elkies' definition -- we require that no quadratic twist
        of a curve L-isogenous to E is a Galois conjugate of E.
        """
        K = E.base_field()
        if K is QQ: return True
        add_field(K, field_label=field_label)

        # first quick test: are the a-invariants Galois invariant?
        if all([is_Galois_invariant(a,field_label) for a in E.ainvs()]):
                return True

        # second quick test: is the conductor invariant?
        if not is_Galois_invariant(E.conductor(),field_label):
                return False

        # third test should catch all non-Q-curves: find primes of
        # good reduction and of the same norm and test if the
        # cardinalities of the reductions are equal.

        pmax = 200
        NN = E.conductor().norm()
        for p in primes(pmax):
            if p.divides(NN):
                continue
            Plist = [P for P in K.primes_above(p)
                     if P.residue_class_degree() == 1]
            if len(Plist)<2:
                continue
            aP = E.reduction(Plist[0]).trace_of_frobenius()
            for P in Plist[1:]:
                if E.reduction(P).trace_of_frobenius() != aP:
                    return False

        # Retrieve precomputed Galois group and isogeny class and test
        G = Glists[K]
        EL = conj_curve(E,G[0]) # base-change to Galois closure
        C = EL.isogeny_class() # cached
        # Here, 'in' does test up to isomorphism!
        return all([conj_curve(E,sigma) in C for sigma in G.gens()])

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
        r"""
        Convert a list of n NF elements to a string consisting of n
        substrings separated by spaces, each substring being a
        comma-separated list of strings representing rational numbers
        representing the NF element with respect to its (power) basis.
        """
        return " ".join([",".join([str(c) for c in list(ai)]) for ai in ainvs])

def ainvs_from_strings(K, ainv_string_list):
        r"""
        Reverse of the previous function: converts a list of strings,
        each representing an NF element, to a list of actual NF
        elements in K.
        """
        return [parse_NFelt(K,ai) for ai in ainv_string_list]

def curve_from_strings(K, ainv_string_list):
        r"""
        Given a number field K and a list of 5 strings, each
        representing an NF element, converts these to elements of K
        and returns the elliptic curve with these a-invariants.
        """
        return EllipticCurve(ainvs_from_strings(K,ainv_string_list))

def minimal_model(E):
        r""" Return a reduced minimal (or semi-minimal) model; here 'reduced'
        means by unit scaling and then by translation.

        NB The isogeny_class function does not currently do any
        minimisation or reduction of models.
        """
        return E.global_minimal_model(E, semi_global=True)

def min_disc_norm(E):
        r"""
        Return the norm of the minimal discriminant ideal of `E`.
        """
        return prod([x.prime().norm()**x.discriminant_valuation() for x in E.local_data()], Integer(1))

# Comparison of curves in one isogeny class using j-invariants, based
# on Lemma: if E1 and E2 are isogenous and not isomorphic (over k)
# then j(E1)!=j(E2) *except* when E1 has potential but not rational CM
# by discriminant d<0 and E2 is the quadratic twist by d of E1.  A
# solution for the tie-break situation was being worked on at ICTP in
# September 2014 by Maarten Derrickx and Heline Deckonick.

# See below this for a key-function version (more efficient, same mathematics)

def curve_cmp(E1,E2):
        r"""
        Comparison function for two isogenous elliptic curves.

        INPUT:

        - ``E1``, ``E2`` -- two isogenous elliptic curves defined over
          the same number field `k`.

        OUTPUT:

        0,+1,-1 (for comparison)
        """
        if E1.is_isomorphic(E2):
                return int(0)

        if E1.has_rational_cm():
                # Order first by the CM discriminant so that curves
                # with the same endo ring are grouped together:
                d1 = E1.cm_discriminant()
                d2 = E2.cm_discriminant()
                t = cmp(d2,d1) # NB the discriminants are negative!
                               # We want -4 before -16
                if t:
                        return t
        # Now either E1 does not have CM or E1, E2 have different endo rings

        # Order by j-invariant if they are different (usually the case):
        c = cmp(E1.j_invariant(),E2.j_invariant())
        if c:
                return c

        # now E1 and E2 must have potential (not rational) CM and be
        # quadratic twists by the CM field.  We do not yet have a good
        # method to distinguish them.  We use the (minimal)
        # discriminant norm, but we do not have a proof that these
        # cannot be equal.
        minD1 = min_disc_norm(E1)
        c = cmp(minD1,min_disc_norm(E2))
        if c:
                return c

        # The above does not distinguish the curves [0,1,0,=3,1] and
        # [0,-i,0,3,i] over Q(i), which are 2-isogenous, i-twists,
        # both have j=8000 and discriminant norm 2^18.  But they can
        # be distinguished by their Tamagawa products (2 and 4),
        # though not by their torsion orders, (both 2).

        tor1 = E1.torsion_order()
        c = cmp(tor1, E2.torsion_order())
        if c:
                return c

        tam1 = E1.tamagawa_product_bsd()
        c = cmp(tam1, E2.tamagawa_product_bsd())
        if c:
                return c

        # todo: resolve the case where E1 and E2 are quadratic twists
        # with the same minimal discriminant norm, torsion order and
        # Tamagawa product!
        print("curves %s and %s are isogenous twists, both with j-invariant %s, minimal discriminant norm %s, torsion order %s, Tamagawa product %s: tie-break condition!" % (E1.ainvs(),E2.ainvs(),E1.j_invariant(),minD1,tor1,tam1))
        return int(0)

# key functions for sorting curves in an isogeny class
def isogeny_class_key_traditional(E):
        return flatten([list(ai) for ai in E.ainvs()])

def isogeny_class_key_cm(E):
        return (int(E.has_rational_cm() and -E.cm_discriminant()),
                flatten([list(ai) for ai in E.ainvs()]))

def isogeny_class_key_j(E):
        return (int(E.has_rational_cm() and -E.cm_discriminant()),
                E.j_invariant(),
                min_disc_norm(E),
                E.torsion_order(),
                E.tamagawa_product_bsd(),
                flatten([list(ai) for ai in E.ainvs()])) # last resort tie-break

isogeny_class_key = isogeny_class_key_j

def Euler_polynomial(E,P):
        r"""
        Return the Euler polynomial of E at the prime P.

        INPUT:

        - `E` -- an elliptic curve defined over a number field K

        - `P` -- a prime ideal of K

        OUTPUT:

        The polynomial `f(X) \in \ZZ[X]` such that `f(N(P)^{-s})` is
        the inverse of the Euler factor of the L-function of `E` at
        `P`.
        """
        EP = E.local_data(P)
        if EP.has_good_reduction():
                return E.reduction(P).frobenius_polynomial().reverse()
        else:
                return 1-EP.bad_reduction_type()*polygen(ZZ)

def rational_Euler_polynomial(E,p):
        r"""
        Return the Euler polynomial of E at the rational prime p.

        INPUT:

        - `E` -- an elliptic curve defined over a number field K

        - `P` -- a prime number

        OUTPUT:

        The polynomial `f(X) \in \ZZ[X]` such that `f(p^{-s})` is
        the inverse of the Euler factor of the L-function of `E` at
        `p`.
        """
        x = polygen(ZZ)
        K = E.base_field()
        return prod([Euler_polynomial(E,P)(x^P.residue_class_degree())
                     for P in K.primes_above(p)])

def rational_L_coefficients(E, nmax, prime_powers_only=True):
        r"""
        Return a dict giving the first ``nmax`` coefficients of the
        L-function of E.

        INPUT:

        - ``E`` -- an elliptic curve defined over a number field.

        - ``nmax`` -- a positive integer

        - ``prime_powers_only`` (bool, default ``True``) -- if
          ``True``, the keys will be restricted to primes powers;
          otherwise all positive integers up to ``nmax``.

        OUTPUT:

        A dict keyed by positive integers `n` up to ``nmax`` whose
        value at `n` is the cofficient of `n^{-s}` in the L-function
        of ``E``, for `n=1,2,\dots,` ``nmax`` (or just the prime
        powers `n>1`).
        """
        # maxexp(p) = max{i: p^i <= nmax}
        lognmax = RR(nmax).log()
        maxexp = lambda p: (lognmax/RR(p).log()).floor()

        polydata = [(p,maxexp(p),rational_Euler_polynomial(E,p))
                    for p in primes(nmax+1)]
        t = PowerSeriesRing(ZZ,'t').gen()
        c = {}
        for p,e,pol in polydata:
                s = (1/pol(t) + O(t**nmax)).dict()
                cp = dict([(p^i,s.get(i,0)) for i in range(1,e+1)])
                c.update(cp)

        # so far, c[n] is defined for n=p^i with 1<n<=nmax, but only when c[n]!=0
        c[1] = Integer(1)
        if prime_powers_only:
                return c

        for n in srange(2,nmax+1):
                if not n in c: # we do not yet know c[n]
                        nf =n.factor()
                        assert len(nf)>1
                        n1 = nf[0][0]**nf[0][1] # p**e
                        n2 = n//n1 # n2<n so we have c[n2]
                        c[n] = c[n1]*c[n2]
        return c

def curve_cmp_via_L(E1,E2,nmax=100):
        r"""
        Comparison function for elliptic curves, using rational L-functions.

        INPUT:

        - ``E1``, ``E2`` - elliptic curves defined over number fields;

        - ``nmax`` (int, default 100) - number of L-series coefficients to use.

        OUTPUT:

        0,+1,-1 (for comparison) based on lexicographical ordering of
        the Dirichlet expansions of the L-functions of E1,E2 (in the
        form `\sum_{n=1}^{\infty}a_n/n^s`).  Since the L-function is
        isogeny-invariant, the output will be 0 only for isogenous
        curves (but see below); this comparison is intended for the
        purpose of sorting isogeny classes.

        .. NOTE:

        If ``nmax`` is too small, the output may be 0 even though the
        curves are not isogenous.
        """
        L1 = rational_L_coefficients(E1,nmax)
        L2 = rational_L_coefficients(E2,nmax)
        L = [L1[n] for n in sorted(L1.keys())]
        c = cmp(L, [L2[n] for n in sorted(L2.keys())])
        if c:
                return c
        # For testing purposes:
        if not E1.is_isogenous(E2):
                print("Warning: curves %s and %s of conductor %s have matching L-functions\n   %s but are not isogenous!" % (E1.ainvs(),E2.ainvs(), E1.conductor(), L))
                return c

# Isogeny class comparison: original form, using the L-functions as
# sums over integral ideals of k.  This matches the sorting of Bianchi
# newforms.

def isog_class_cmp1(k, I, J):
	E_I = curve_from_strings(k,I[0].split()[6:11])
	E_J = curve_from_strings(k,J[0].split()[6:11])

	for p in Plists[k]:
		c = int(ap(E_I, p) - ap(E_J, p))
		if c: return cmp(c,0)

	raise NotImplementedError("Bound on primes is too small to determine...")

# Isogeny class comparison: experimental for, based on comparison
# between the L-functions as rational Dirichlet series (indexed by
# postive integers, hence no choices needed!)  implemented at ICTP
# September 2014 by Angelos Koutsianas, Alejandro Argaez, Daniel
# Kohen, and revised here by John Cremona.  NB This has been
# *discarded* since Galois conjugate classes have the same rational
# L-function and would need a tie-break anyway.

def isog_class_cmp2(k, I, J):
	E1 = curve_from_strings(k,I[0].split()[6:11])
	E2 = curve_from_strings(k,J[0].split()[6:11])
        return curve_cmp_via_L(E1,E2)

fields = {} # keys are field labels, values are NumberFields
import yaml
field_dict = yaml.load(file("HMField_data.yaml")) # all the totally real fields in the LMFDB

def field_from_label(lab):
        if lab in fields:
                return fields[lab]
        dummy, deg, sig, abs_disc = field_data(lab)
        x = polygen(QQ)
        name = 'a'
        if deg==2:
                d = ZZ(abs_disc)
                if d==4: name='i'
                if d==8: name='t'
                if d==3: name='w'
                if sig[0]==0: d=-d
                t = d%4
                assert t in [0,1]
                pol = x^2 - t*x + (t-d)/4
        elif lab=='3.1.23.1':
                pol = x**3 - x**2 +1
        else:
            if lab in field_dict:
                coeffs = field_dict[lab]['coeffs']
                pol = x.parent()(coeffs)
            else:
                raise NotImplementedError("cannot yet handle field %s" % lab)
        K = NumberField(pol, name)
        fields[lab] = K
        print "Created field from label %s: %s" % (lab,K)
        return K

def read_curves(infile, only_one=False, ncurves=0):
        r"""
        Iterator to loop through lines of a curves.* file each
        containing 13 data fields as defined the in the ecnf-format.txt file,
        yielding its curves as EllipticCurve objects.

        If only_one is True, skips curves whose 4th data field is
        *not* 1, hence only yielding one curve per isogeny class.
        """
        count=0
        for L in file(infile).readlines():
                #sys.stdout.write(L)
                data = L.split()
                if len(data)!=13:
                        print "line %s does not have 13 fields, skipping" % L
                        continue
                if only_one and data[3]!='1':
                        continue
                count +=1
                if ncurves and count>ncurves:
                       raise StopIteration
                field_label = data[0]
                K = field_from_label(field_label)
                N_label = data[1]
                iso_label = data[2]
                c_num = data[3]
                N_def = data[4]
                E = curve_from_strings(K, data[6:11])
                yield (field_label,N_label,N_def,iso_label,c_num,E)

def get_isoclass_letter(N, E):
        K = E.base_field()
        aplist = [ap(E,P) for P in Plists[K]][:25]
        #print("N=%s, ap=%s" % (N,aplist))
        try:
                nfs = nf_data[K][N]
        except KeyError: # No newform of this level at all
                #print("(1) No newform matches %s, conductor %s" % (E.ainvs(),N))
                n = cm_counts[K].get(N,0)
                cm_counts[K][N] = n+1
                return "CM"+cremona_letter_code(n)
        #print("Trying newforms %s" % nfs.keys())
        for id in nfs.keys():
                nfap = nfs[id]['ap'][:25]
                if nfap == aplist:
                        #print("Matches %s" % id)
                        return id
                #else: print("No match with %s: ap=%s" % (id,nfap))

        # If we get to here, no newform of this level matches the ap
        #print("(2) No newform matches %s, conductor %s" % (E.ainvs(),N))
        n = cm_counts[K].get(N,0)
        cm_counts[K][N] = n+1
        return "CM"+cremona_letter_code(n) # 0->a, 1->b etc

# Main curve processing function
def process_curves(curves, outfile = None, classfile=None, verbose=0):
        r"""
        Given a list or iterator yielding a sequence of elliptic
        curves (which could be either an actual list of elliptic
        curves, or something like read_curves(file_name)), processses
        these and writes the results to an output file (if given) and
        to the screen (if verbose>0).

        The input curves do not have to be defined over the same base
        field; the output will be sorted first by field.

        Each curve is first compared with a list of curves previously
        encountered to see if it is isomorphic *or isogenous* to any
        of these, in which case it is ignored.  Hence, if the input
        file contains several curves in an isogeny class, all but the
        first will effectively be ignored.  After that the complete
        isogeny class is computed, sorted, and data for output
        computed for each curve.

        Finally, after all the input and processing are complete, the
        whole lot is output to the file and/or screen, sorted as
        follows: by field, then conductor norm, then conductor (sorted
        using the HNF of the ideal), then by isogeny class, then by
        curves in the class.  By default, the letter labels for
        isogeny classes are generated on the fly, which should be fine
        provided that the input is complete (at least one curve in
        every isogeny class) for each conductor present.  Optionally,
        they can be generated by comparison with the associated
        newform: currently only available for Bianchi newforms over
        K=Q(sqrt(-d)) for d=1,2,3,7,11.  In these cases, since curves
        with CM rational over K are not attached to Bianchi newforms,
        the letter labels are of the form "CMa", "CMb", etc and not
        just "a", "b", etc.
        """
        if outfile:
                outfile = file(outfile, mode="a")
        if classfile:
                classfile = file(classfile, mode="a")

	data = {} # holds line to be output into the main output file,
                  # indexed by (1) field (2) conductor norm (3,4) HNF
                  # generators (only sufficient for quadratic fields!)
                  # i.e. data[k][n][c][d] is a list of strings to be
                  # output, each defining one curve over field k,
                  # conductor norm n, conductor HNF <n/d,c+d*alpha>
                  # where Z[alpha] is the ring of integers.

        isogdata = {} # ditto for isogeny classes & isog matrix data file.

	for E in curves:
                N_label = None
                field_label = None
                if isinstance(E,tuple):
                        field_label, N_label, N_def, iso_label, c_num, E = E
                        if verbose>0:
                                print("processing E = %s with conductor %s, label %s-%s%s... over field %s" % (list(E.ainvs()),N_def,N_label,iso_label,c_num,field_label))
                else:
                        if verbose>0:
                                print("processing E = %s..." % list(E.ainvs()))
                k = E.base_field()
                add_field(k, field_label=field_label)
                D = Dlists[k]
                G = Glists[k]
                used = used_curves[k]
                isog_class_cmp = ic_cmp[k]
                field_label = labels[k]
                nf_data_k = nf_data[k] # may be None; that's ok
                # Set up 2 dicts to collect the data to be output after sorting:
                if not k in data:
                        data[k] = {}
                        isogdata[k] = {}
                data_k = data[k]
                isogdata_k = isogdata[k]

                E = minimal_model(E)
		N = E.conductor()
		norm = N.norm()
                if N_label:
                    if not norm==ZZ(N_label.split(".")[0]):
                        print("****************************")
                        print("Conductor norm = %s, does not match label %s" % (norm, N_label))
                        print("****************************")

                E2 = found(E, norm)
		if E2:
                    if verbose>0:
                        print(" -- isogenous to previous curve %s" % list(E2.ainvs()))
                else:
                        if verbose>1:
                                print(" -- new isogeny class")
			# Conductor
			hnf = N.pari_hnf()
                        if N_label:
                            cond_label = N_label
                            cond_def = N_def
                        else:
                            cond_def = "[%i,%s,%s]" % (norm, hnf[1][0], hnf[1][1])
                            cond_label = cond_def
                        if nf_data_k:
                            isog = get_isoclass_letter(cond_label,E)
                            if verbose>0:
                                print("Isogeny class label = %s.%s" % (cond_label,isog))
                        else:  # default
                            if N_label:
                                isog = iso_label
                            else:
                                isog = ":isog"

			# Setup data
			if norm not in data_k:
				data_k[norm] = {}
				isogdata_k[norm] = {}
			if hnf[1][0] not in data_k[norm]:
				data_k[norm][hnf[1][0]] = {}
				isogdata_k[norm][hnf[1][0]] = {}
			if hnf[1][1] not in data_k[norm][hnf[1][0]]:
				data_k[norm][hnf[1][0]][hnf[1][1]] = []
				isogdata_k[norm][hnf[1][0]][hnf[1][1]] = []
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
					if E.is_isogenous(curve_from_strings(k, curve_data[6:11]), maxnorm=200):
                                                print("Warning: input curve %s isogenous to previous curve %s but not found by isogeny class computation!" % (list(E.ainvs()),curve_data[6:11]))
						curve_data[3] = len(found_isog_class)+1
						curve_data[6:11] = [",".join([str(c) for c in ai]) for ai in ainvs]
						data_k[norm][hnf[1][0]][hnf[1][1]][n].append(" ".join(curve_data))
						break

			# Compute the isogeny class
                        if verbose>1:
                                print("computing the isogeny class")
			Cl = E.isogeny_class()
                        clist0 = [minimal_model(C) for C in Cl.curves]
                        mat0 = Cl.matrix()
                        # sort into new order (will be redundant later)
                        clist = sorted(clist0, key=isogeny_class_key)
                        # perm[i]=j where sorted#i = unsorted#j
                        perm = dict([(i,clist0.index(E)) for i,E in enumerate(clist)])
                        mat = copy(mat0) # to set the size etc
                        for i in range(len(clist)):
                                for j in range(len(clist)):
                                        mat[i,j] = mat0[perm[i],perm[j]]

			if norm not in used:
				used[norm] = []
                        if verbose and len(clist)>1:
                            print("Adding %s isogenous curves" % (len(clist)-1))
                            #for E2 in clist[1:]:
                            #    print list(E2.ainvs())
			used[norm] += clist

                        matlist = str([list(r) for r in mat]).replace(' ','')
                        isogdata_line = "%s %s %s %i %s" % (field_label, cond_label, isog, 1, matlist)
			isogdata_k[norm][hnf[1][0]][hnf[1][1]].append([isogdata_line])
                        if verbose>1:
                                print("isog_data: %s" % isogdata_line)

                        # Q-curve? (isogeny class invariant)
                        q_curve = int(is_Q_curve(E))
                        if verbose>1:
                            if q_curve:
                                print("...Q-curve")
                            else:
                                print("...not a Q-curve")

			tmp = [] # list of output lines (with
                                 # placeholder for isog code, filled
                                 # in after sorting)

			for n, E2 in enumerate(clist):
				# a-invs
				ainvs = E2.a_invariants()
                                ainv_string = ainvs_to_string(ainvs)
				# Disc
				j = E2.j_invariant()
				disc = cm_j_invariants.get(j, 0)

                                curve_line = "%s %s %s %i %s %i %s %i %i" % (field_label, cond_label, isog, n + 1, cond_def, norm, ainv_string, disc, q_curve)
                                if verbose>1:
                                        print("curve_line: %s" % curve_line)
				tmp.append(curve_line)
			data_k[norm][hnf[1][0]][hnf[1][1]].append(tmp)

	# Sort and output the data

        ks = data.keys()
        if verbose>0:
                print
                print "fields: %s" % ks
        ks.sort()
        for k in ks:
            data_k = data[k]
            isogdata_k = isogdata[k]
            norms = data_k.keys()
            norms.sort()
            for norm in norms:
                data_k_n = data_k[norm]
                isogdata_k_n = isogdata_k[norm]
		hnf0s = data_k_n.keys()
		hnf0s.sort()
		for hnf0 in hnf0s:
                        data_k_n_h = data_k_n[hnf0]
                        isogdata_k_n_h = isogdata_k_n[hnf0]
			hnf1s = data_k_n_h.keys()
			hnf1s.sort()
			for hnf1 in hnf1s:
                                dat = data_k_n_h[hnf1]
                                isogdat = isogdata_k_n_h[hnf1]
				dat.sort(cmp = isog_class_cmp)
				for n, (cdata,isodata) in enumerate(zip(dat,isogdat)):
                                        isoline = isodata[0]
                                        if isog==":isog":
                                                isog_letter = cremona_letter_code(n)
                                                isoline = isoline.replace(":isog", isog_letter)
                                        if classfile:
                                                classfile.write(isoline+'\n')
                                        if verbose>0:
                                                print isoline
					for E_data in cdata:
                                                line = E_data
                                                if isog==":isog":
                                                        isog_letter = cremona_letter_code(n)
                                                        line = line.replace(":isog", isog_letter)
                                                if outfile:
                                                        outfile.write(line+'\n')
						if verbose>0:
                                                        print line

# Old version of previous function, will become redundant.

# Basic info about the curves
def basic_info(curves, outfile = None, classfile=None, verbose=0):
        r"""
        Given a list or iterator yielding a sequence of elliptic
        curves (which could be either an actual list of elliptic
        curves, or something like read_curves(file_name)), processses
        these and writes the results to an output file (if given) and
        to the screen (if verbose>0).

        The input curves do not have to be defined over the same base
        field; the output will be sorted first by field.

        Each curve is first compared with a list of curves previously
        encountered to see if it is isomorphic *or isogenous* to any
        of these, in which case it is ignored.  Hence, if the input
        file contains several curves in an isogeny class, all but the
        first will effectively be ignored.  After that the complete
        isogeny class is computed, sorted, and data for output
        computed for each curve.

        Finally, after all the input and processing are complete, the
        whole lot is output to the file and/or screen, sorted as
        follows: by field, then conductor norm, then conductor (sorted
        using the HNF of the ideal), then by isogeny class (with
        letter labels created on the fly after sorting), then by
        curves in the class.
        """
        if outfile:
                outfile = file(outfile, mode="a")
        if classfile:
                classfile = file(classfile, mode="a")

	data = {} # holds line to be output into the main output file,
                  # indexed by (1) field (2) conductor norm (3,4) HNF
                  # generators (only sufficient for quadratic fields!)
                  # i.e. data[k][n][c][d] is a list of strings to be
                  # output, each defining one curve over field k,
                  # conductor norm n, conductor HNF <n/d,c+d*alpha>
                  # where Z[alpha] is the ring of integers.

        isogdata = {} # ditto for isogeny classes & isog matrix data file.

	for E in curves:
                N_label = None
                field_label = None
                if isinstance(E,tuple):
                        N_label, E = E
                        if verbose>0:
                                print("processing E = %s with conductor %s..." % (list(E.ainvs()),N_label))
                else:
                        if verbose>0:
                                print("processing E = %s..." % list(E.ainvs()))
                k = E.base_field()
                add_field(k, field_label=field_label)
                D = Dlists[k]
                G = Glists[k]
                used = used_curves[k]
                isog_class_cmp = ic_cmp[k]
                field_label = labels[k]
                if not k in data:
                        data[k] = {}
                        isogdata[k] = {}
                data_k = data[k]
                isogdata_k = isogdata[k]

		# Get a global minimal model for E if possible
                E = E.global_minimal_model(semi_global=True)
		N = E.conductor()
		norm = N.norm()

		if found(E, norm):
                        if verbose>0:
                                print(" -- isogenous to a previous curve")
                else:
                        if verbose>1:
                                print(" -- new isogeny class")
			# Conductor
			hnf = N.pari_hnf()
                        if N_label:
                                cond_label = N_label
                        else:
                                cond_label = "[%i,%s,%s]" % (norm, hnf[1][0], hnf[1][1])

			# Setup data
			if norm not in data_k:
				data_k[norm] = {}
				isogdata_k[norm] = {}
			if hnf[1][0] not in data_k[norm]:
				data_k[norm][hnf[1][0]] = {}
				isogdata_k[norm][hnf[1][0]] = {}
			if hnf[1][1] not in data_k[norm][hnf[1][0]]:
				data_k[norm][hnf[1][0]][hnf[1][1]] = []
				isogdata_k[norm][hnf[1][0]][hnf[1][1]] = []
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
                                                print("Warning: input curve %s isogenous to a previous curve but not found by isogeny class computation!" % list(E.ainvs()))
						curve_data[3] = len(found_isog_class)+1
						curve_data[6:11] = [",".join([str(c) for c in ai]) for ai in ainvs]
						data_k[norm][hnf[1][0]][hnf[1][1]][n].append(" ".join(curve_data))
						break

			# Find the isogeny class
                        if verbose>1:
                                print("computing the isogeny class")
			Cl = E.isogeny_class()
                        clist0 = [minimal_model(C) for C in Cl.curves]
                        mat0 = Cl.matrix()
                        # sort into new order (will be redundant later)
                        clist = sorted(clist0, cmp=curve_cmp)
                        perm = dict([(i,clist0.index(E)) for i,E in enumerate(clist)])
                        mat = copy(mat0) # to set the size etc
                        for i in range(len(clist)):
                                for j in range(len(clist)):
                                        mat[perm[i],perm[j]] = mat0[i,j]

			if norm not in used:
				used[norm] = []
			used[norm] += clist

                        matlist = str([list(r) for r in mat]).replace(' ','')
                        isogdata_line = "%s %s :isog %i %s" % (field_label, cond_label, 1, matlist)
			isogdata_k[norm][hnf[1][0]][hnf[1][1]].append([isogdata_line])
                        if verbose>1:
                                print("%s" % isogdata_line)

                        # Q-curve? (isogeny class invariant)
                        q_curve = int(is_Q_curve(E))

			tmp = [] # list of output lines (with
                                 # placeholder for isog code, filled
                                 # in after sorting)

			for n, E2 in enumerate(clist):
				# a-invs
				ainvs = E2.a_invariants()
                                ainv_string = ainvs_to_string(ainvs)
				# Disc
				j = E2.j_invariant()
				disc = cm_j_invariants.get(j, 0)

				tmp.append("%s %s :isog %i %s %i %s %i %i" % (field_label, cond_label, n + 1, cond_label, norm, ainv_string, disc, q_curve))
                        #print "appending %s curves" % len(tmp)
			data_k[norm][hnf[1][0]][hnf[1][1]].append(tmp)

	# Sort and output the data

        ks = data.keys()
        if verbose>0:
                print
                print "fields: %s" % ks
        ks.sort()
        for k in ks:
            data_k = data[k]
            isogdata_k = isogdata[k]
            norms = data_k.keys()
            norms.sort()
            for norm in norms:
                data_k_n = data_k[norm]
                isogdata_k_n = isogdata_k[norm]
		hnf0s = data_k_n.keys()
		hnf0s.sort()
		for hnf0 in hnf0s:
                        data_k_n_h = data_k_n[hnf0]
                        isogdata_k_n_h = isogdata_k_n[hnf0]
			hnf1s = data_k_n_h.keys()
			hnf1s.sort()
			for hnf1 in hnf1s:
                                dat = data_k_n_h[hnf1]
                                isogdat = isogdata_k_n_h[hnf1]
				dat.sort(cmp = isog_class_cmp)
				for n, (cdata,isodata) in enumerate(zip(dat,isogdat)):
					isog_letter = cremona_letter_code(n)
                                        isoline = isodata[0].replace(":isog", isog_letter)
                                        if classfile:
                                                classfile.write(isoline+'\n')
                                        if verbose>0:
                                                print isoline
					for E_data in cdata:
                                                line = E_data.replace(":isog", isog_letter)
                                                if outfile:
                                                        outfile.write(line+'\n')
						if verbose>0:
                                                        print line

def run1(pth,fld):
    infile = "%s/curves1.%s" % (pth,fld)
    outfile = "%s/curves.%s" % (pth,fld)
    classfile = "%s/isoclass.%s" % (pth,fld)
    process_curves(read_curves(infile), outfile=outfile, classfile=classfile, verbose=1)

def run2(pth,fld):
    infile = "%s/curves2.%s" % (pth,fld)
    outfile = "%s/curves.%s" % (pth,fld)
    classfile = "%s/isoclass.%s" % (pth,fld)
    process_curves(read_curves(infile), outfile=outfile, classfile=classfile, verbose=1)
