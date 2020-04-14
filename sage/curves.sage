# coding=utf-8

from sage.all import cm_j_invariants_and_orders, ZZ, QQ, RR, EllipticCurve, flatten, legendre_symbol, polygen, prod, primes, PowerSeriesRing, Integer, NumberField, srange, copy, IntegerRing, Primes, O, Magma
from sage.databases.cremona import cremona_letter_code
from psort import (nf_key, primes_of_degree_iter, primes_iter, ideal_label)
from nfscripts import ideal_HNF
from fields import nf_lookup

# Originally adapted from Warren Moore's scripts.  The main function
# process_curves(curves) takes a list of elliptic curves and outputs
# data files.  Isogenous curves are computed and sorted.

# NB requires functionality of Sage-7.1.

# cached field data: these will all be keyed by fields k, values as shown:
Plists = {} # sorted list of primes of k
Dlists = {} # |disc(k)|
Glists = {} # Gal(k/Q)
autolists = {} # automorphisms of k
galois_flags = {} # is k Galois?
labels = {} # label of k
pols   = {} # defining polynomial of k as comma-sep. string
nf_data = {} # newform data for k (if available)
used_curves = {} # dict with key by norm(conductor), value a list of
                 # curves so far processed witth that norm-conductor
class_key = {} # dict whose values are the isogeny class sort key
               # function (to be used with curves of the same
               # conductor only)
cm_counts = {} # dict with keys conductor labels, values counts of
               # classes with rational CM (imaginary quadratic fields
               # only) for labelling of these, as they do not have
               # associated Bianchi newforms
#
cm_j_dict = {} # keys are CM j-invariants, values are associated discriminanrs


def add_field(K, field_label=None, prime_norm_bound=200):
    if K in used_curves:
        return

    print("Adding {} (label={}) to fields collection...".format(K,field_label))
    Plists[K] = list(primes_iter(K,maxnorm=prime_norm_bound))
    absD = K.discriminant().abs()
    s = K.signature()[0] # number of real places
    d = K.degree()
    pols[K] = ",".join([str(c) for c in list(K.defining_polynomial())])

# Warning: number fields whose label's 4'th component is not 1 will
# not be handled correctly here; this is only an issue when there is
# more than one field of given signature and abs(disc), so fine for
# quadratic fields.  This problem first hit for field 4.4.16448.2.
# When we processed the curves for that field they were in a different
# input file so were not mixed up with curves for 4.4.16448.1 luckily,
# and manual editing of the output was sufficient.

    if field_label==None:
        field_label = '%s.%s.%s.1' % (str(d),str(s),str(absD))
        print("...created new label {}".format(field_label))
    labels[K] = field_label
    Dlists[K] = absD
    autolists[K] = K.automorphisms()
    galois_flags[K] = len(autolists[K])==K.degree()
    if d<5:
        Glists[K] = K.galois_group(names='b')
    #print("...finding and caching Galois closure...")
    #K.SF = K.defining_polynomial().splitting_field('b')
    #print("...finding CM j-invariants...")
    for dd, f, j in cm_j_invariants_and_orders(K):
	cm_j_dict[j] = dd * (f ^ 2)
    used_curves[K] = {}
    nf_data[K] = None
    cm_counts[K] = {}
    if d==2 and s==0 and absD in [3,4,7,8,11]:
            from nfscripts import read_newform_data
            #from nfscripts import nf_filename_from_D
            #nf_filename = nf_filename_from_D(absD)
            #nf_filename = "/home/jec/bianchi-data/nflist/nflist.11.20001-30000"
            #nf_filename = "/home/jec/bianchi-data/newforms/newforms.3.120001-130000"
            nf_filename = "/home/jec/bianchi-data/newforms/newforms.3.missing_curves"
            print("reading newform data from {}".format(nf_filename))
            nf_data[K] = read_newform_data(nf_filename)
    print("...finished adding field.")

# utiity functions for converting ideal labels to LMFDB standard

the_ideal_labels = {}

def convert_ideal_label(K, lab, IQF_only=True):
    """An ideal label of the form N.c.d is converted to N.i.  Here N.c.d
    defines the ideal I with Z-basis [a, c+d*w] where w is the standard
    generator of K, N=N(I) and a=N/d.  The standard label is N.i where I is the i'th ideal of norm N in the standard ordering.

    NB Only intended for use in coverting IQF labels!  To get the standard label from any ideal I just use ideal_label(I).

    """
    global the_ideal_labels
    if K in the_ideal_labels:
        if lab in the_ideal_labels[K]:
            return the_ideal_labels[K][lab]
        else:
            pass
    else:
        the_ideal_labels[K] = {}

    comps = lab.split(".")
    # test for labels which do not need any conversion
    if len(comps)==2:
        return lab
    assert len(comps)==3
    N, c, d = [int(x) for x in comps]
    a = N//d
    I = K.ideal(a, c+d*K.gen())
    newlab = ideal_label(I)
    #print("Ideal label converted from {} to {} over {}".format(lab,newlab,K))
    the_ideal_labels[K][lab] = newlab
    return newlab


def convert_conductor_label(field_label, label):
    """If the field is imaginary quadratic, calls convert_ideal_label, otherwise just return label unchanged.
    """
    if field_label.split(".")[:2] != ['2','0']:
        return label
    K = nf_lookup(field_label)
    new_label = convert_ideal_label(K,label)
    #print("Converting conductor label from {} to {}".format(label, new_label))
    return new_label

def ap(E, p):
        r"""
        Return a_p(E).

        INPUT:

        - ``E`` - an elliptic curve defined over a number field `k`;

        - ``p`` - a prime ideal of `k`.

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

def ap_list(E):
        r"""
        Return [a_p(E) for p in Plists[k]] where k=E.base_field().

        INPUT:

        - ``E`` - an elliptic curve defined over a number field `k`;

        OUTPUT:

        A list of a_P(E) for P in the standard, ordered list
        Plists[k], where k is the base field of E.
        """
        K = E.base_field()
        add_field(K)
        return [ap(E,p) for p in Plists[K]]

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

def is_Q_curve(E, field_label=None, verbose=False):
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
        if verbose:
            print("Checking whether {} is a Q-curve".format(E))

        # all curves with rational j-invariant are Q-curves:
        jE = E.j_invariant()
        if jE in QQ:
            if verbose:
                print("Yes: j(E) is in QQ")
            return True

        K = E.base_field()
        jpoly = jE.minpoly()
        if jpoly.degree()<K.degree():
            print("switching to smaller base field: j's minpoly is {}".format(jpoly))
            K1, _ = K.subfield(jE, 'j')
            #print("K1 = {}".format(K1))
            K2, iso, inv = K1.optimized_representation()
            #print("K2 = {}".format(K2))
            jE = inv(K1.gen())
            print("New j is {} with minpoly {}".format(jE, jE.minpoly()))
            assert jE.minpoly()==jpoly
            E = EllipticCurve(j=jE)
            K = K2
            field_label=None
            print("New test curve is {}".format(E))

        # CM curves are Q-curves:
        if E.has_cm():
            if verbose:
                print("Yes: E is CM")
            return True

        add_field(K, field_label=field_label)

        # Simple test should catch many non-Q-curves: find primes of
        # good reduction and of the same norm and test if the
        # traces of Frobenius are equal *up to sign*

        pmax = 1000
        NN = E.conductor().norm()
        for p in primes(pmax):
            if p.divides(NN):
                continue
            Plist = [P for P in K.primes_above(p)
                     if P.residue_class_degree() == 1]
            if len(Plist)<2:
                continue
            aP0 = E.reduction(Plist[0]).trace_of_frobenius()
            for P in Plist[1:]:
                aP = E.reduction(P).trace_of_frobenius()
                if aP.abs() != aP0.abs():
                    if verbose:
                        print("No: incompatible traces of Frobenius at primes above {}: {} and {}".format(p,aP0,aP))
                    return False

        if verbose:
            print("...all aP test pass for p<{}".format(pmax))

        C = E.isogeny_class()
        jC = [E2.j_invariant() for E2 in C]
        if any(j in QQ for j in jC):
            if verbose:
                print("j-invariants in class: {}".format(jC))
                print("Yes: an isogenous curve has j in QQ")
            return True

        # Galois case:
        if galois_flags[K]:
            autos = autolists[K]
            t = all(sigma(jE) in jC for sigma in autos)
            if verbose:
                print("j-invariants in class: {}".format(jC))
                if t:
                    print("Yes: class contains all conjugate j-invariants")
                else:
                    print("No: class does not contain all conjugate j-invariants")
            return t

        if verbose:
            print("K not Galois, E might be a Q-curve but not yet proved")

        if K.degree() in [3,4]:
            L = K.galois_closure('b')
            emb = K.embeddings(L)[0]
            EL = E.change_ring(emb)
            if verbose:
                print("Base changing to the Galois closure {}".format(L))
            # Don't do this as it will cause an infinite loop
            # return is_Q_curve(EL, verbose=verbose)
            autos = L.automorphisms()
            jE = emb(jE)
            nC = len(C)
            from sage.schemes.elliptic_curves.gal_reps_number_field import reducible_primes_naive
            red_pr = reducible_primes_naive(EL, 31, 100)
            # Compute partial isogeny class only using p-iogenies
            # for the primes p in red_pr.  This is much faster since
            # the slow part of computing isgeny classes is finding
            # the reducible primes.
            if verbose:
                print("Computing isogeny class using only the primes in {}".format(red_pr))
            CL = EL.isogeny_class(red_pr)
            if len(CL)==nC:
                if verbose:
                    print(" -- no more curves in the class")
            else:
                nC=len(CL)
                if verbose:
                    print(" -- class now contains {} curves".format(nC))
                jCL = [E2.j_invariant() for E2 in CL]
                t = all(sigma(jE) in jCL for sigma in autos)
                if verbose:
                    print("j-invariants in class: {}".format(jCL))
                    if t:
                        print("Yes: class contains all conjugate j-invariants")
                    else:
                        print("No: class does not contain all conjugate j-invariants")
                if t:
                    return t

            if verbose:
                print("%%%%%%%%%%")
                print("No conclusion after computing the isogeny class over the Galois closure using these primes")
                print("%%%%%%%%%%")
                print("Computing full isogeny class...")
            CL = EL.isogeny_class()
            if len(CL)==nC:
                if verbose:
                    print("No more isogenies, so not a Q-curve")
                return False
            print("Full isogeny class has {} curves".format(len(CL)))
            jCL = [E2.j_invariant() for E2 in CL]
            t = all(sigma(jE) in jC for sigma in autos)
            if verbose:
                print("j-invariants in class: {}".format(jC))
                if t:
                    print("Yes: class (using first {} primes) contains all conjugate j-invariants".format(n+1))
                else:
                    print("No: class (using first {} primes) does not contain all conjugate j-invariants".format(n+1))
            return t

        return '?'

        # # Retrieve precomputed Galois group and isogeny class and test
        # G = Glists[K]
        # EL = conj_curve(E,G[0]) # base-change to Galois closure
        # C = EL.isogeny_class() # cached
        # # Here, 'in' does test up to isomorphism!
        # return all([conj_curve(E,sigma) in C for sigma in G.gens()])

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
        I = E.minimal_discriminant_ideal()
        if I.ring()==IntegerRing():
            return I.gen()
        return I.norm()

# Comparison of curves in one isogeny class using j-invariants, based
# on Lemma: if E1 and E2 are isogenous and not isomorphic (over k)
# then j(E1)!=j(E2) *except* when E1 has potential but not rational CM
# by discriminant d<0 and E2 is the quadratic twist by d of E1.  One
# solution for the tie-break situation was being worked on at ICTP in
# September 2014 by Maarten Derrickx and Heline Deckonick.  The
# solution implemented here was developed by Andrew Sutherland and
# John Cremona in June 2016, and requires using the "first" degree 1
# prime with certain properties, hence requires a fixed ordering of
# prime ideals.  This was also developed by Andrew Sutherland, John
# Cremona and Aurel Page in June 2016.


# key functions for sorting curves in an isogeny class
def isogeny_class_key_traditional(E):
        return flatten([list(ai) for ai in E.ainvs()])

def isogeny_class_key_cm(E):
        return (int(E.has_rational_cm() and -E.cm_discriminant()),
                flatten([list(ai) for ai in E.ainvs()]))

# A version of primes_of_degree_iter for K=Q:
def primes_iter_Q(condition):
    for p in Primes():
        if condition(p):
            yield(p)

def cmj_key(E):
    r""" Key to compare curves with non-rational CM which are quadratic
    twists over the CM field.  This will be called on lots of curves
    for which this tie-break comparison is not needed, so we return 0
    instantly when we know that is the case.
    """
    if (not E.has_cm()) or E.has_rational_cm():
        return 0
    d = E.cm_discriminant()
    K = E.base_field()
    deg = K.absolute_degree()
    D = 1 if deg==1 else ZZ(K.defining_polynomial().discriminant())
    j = E.j_invariant()
    c4, c6 = E.c_invariants()
    jj, c, w = (j, c4, 4) if j==1728 else (j-1728, c6, 6)
    NN = E.conductor() if deg==1 else E.conductor().norm()
    bad = 6*d*D*NN

    # Get the first degree 1 prime P, dividing a prime p not dividng
    # bad for which d is a quadratic non-residue, such that j-1728 (or
    # j when j=1728) is a P-unit:
    ptest = lambda p: not p.divides(bad) and legendre_symbol(d,p)==-1
    if deg==1:
        it = primes_iter_Q(ptest)
    else:
        it = primes_of_degree_iter(K,deg=1, condition = ptest)
    P = it.next()
    while jj.valuation(P)!=0:
        P = it.next()
    p = P if deg==1 else P.smallest_integer() # = residue characteristic
    print("E = {} with j = {}: tie-break prime P = {} above p = {}".format(E.ainvs(), j, P, p))

    # The key is now (c6|p) (or (c4|p) if c6=0) with c4, c6 from the
    # P-minimal model.  Although E has good reduction at P the model
    # may not be minimal, and some adjustment is necessary:
    k = c.valuation(P)
    if k>0:
        assert w.divides(k)
        pi = K.uniformizer(P,others='negative')
        c = c/pi**(k//w) # still integral everywhere
    return legendre_symbol(K.residue_field(P)(c),p)

def isomorphism_class_key_j(E):
    """FOr isogenous curves, first sort by CM-discriminant, then by
    j-invariant, then (only necessary when E has potential CM) the
    tie-break.

    """
    return (int(E.has_rational_cm() and -E.cm_discriminant()),
            nf_key(E.j_invariant()),
            cmj_key(E))

isomorphism_class_key = isomorphism_class_key_j

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

def read_curve_file(infile):
    r"""
    Read a curves file, each line containing 13 data fields as defined
    the in the ecnf-format.txt file. Output is a list of dicts, one
    per curve, in the same order as input.
    """
    curves = []
    index = 0
    for L in file(infile).readlines():
        if L[0]=='#': # allow for comment lines
            continue
        data = L.split()
        if len(data)!=13:
            print("line {} does not have 13 fields, skipping".format(L))
            continue
        index += 1
        curve = {'index': index,
                 'field_label': data[0],
                 'N_label': data[1],
                 'iso_label': data[2],
                 'c_num': data[3],
                 'N_def': data[4],
                 'N_norm': data[5],
                 'ainvs': data[6:11],
                 'cm_flag': data[11],
                 'q_curve_flag': data[12]
                 }
        curves.append(curve)
    curves.sort(key = lambda c: c['index'])
    return curves

def read_curvedata_file(infile):
    r"""
    Read a curvedata file, each line containing 9+ data fields as defined
    the in the ecnf-format.txt file. Output is a list of dicts, one
    per curve, in the same order as input.
    """
    curves = []
    index = 0
    for L in file(infile).readlines():
        if L[0]=='#': # allow for comment lines
            continue
        data = L.split()
        ngens = int(data[7])
        if len(data) != 9+ngens:
            print("line {} does not have 9+{} fields, skipping".format(L, ngens))
            continue
        index += 1
        curve = {'index': index,
                 'field_label': data[0],
                 'N_label': data[1],
                 'iso_label': data[2],
                 'c_num': data[3],
                 'rank': data[4],
                 'rank_bounds': data[5],
                 'an_rank': data[6],
                 'ngens': data[7],
                 'gens': data[8:8+ngens],
                 'sha_an': data[8+ngens]
                 }
        curves.append(curve)
    curves.sort(key = lambda c: c['index'])
    return curves

def read_isoclass_file(infile):
    r"""
    Read an isoclass file, each line containing 5 data fields as defined
    the in the ecnf-format.txt file. Output is a list of dicts, one
    per curve, in the same order as input.
    """
    curves = []
    index = 0
    for L in file(infile).readlines():
        if L[0]=='#': # allow for comment lines
            continue
        data = L.split()
        if len(data) != 5:
            print("line {} does not have 5 fields, skipping".format(L))
            continue
        index += 1
        curve = {'index': index,
                 'field_label': data[0],
                 'N_label': data[1],
                 'iso_label': data[2],
                 'c_num': data[3],
                 'isomat': data[4]
                 }
        curves.append(curve)
    curves.sort(key = lambda c: c['index'])
    return curves

def write_curve_file(curves, outfile):
    r"""
    Write a curves file, each line containing 13 data fields as defined
    the in the ecnf-format.txt file. Input is a list of dicts.
    """
    out = file(outfile, 'w')
    for c in curves:
        line = " ".join([c['field_label'],
                         c['N_label'],
                         c['iso_label'],
                         c['c_num'],
                         c['N_def'],
                         c['N_norm']] + c['ainvs'] + [c['cm_flag'],
                         c['q_curve_flag']])
        out.write(line+"\n")
    out.close()

def write_curvedata_file(curves, outfile):
    r"""
    Write a curves file, each line containing 9+ngens data fields as defined
    the in the ecnf-format.txt file. Input is a list of dicts.
    """
    out = file(outfile, 'w')
    for c in curves:
        line = " ".join([c['field_label'],
                         c['N_label'],
                         c['iso_label'],
                         c['c_num'],
                         c['rank'],
                         c['rank_bounds'],
                         c['an_rank'],
                         c['ngens']] + c['gens'] + [c['sha_an']])
        out.write(line+"\n")
    out.close()

def write_isoclass_file(curves, outfile):
    r"""
    Write an isoclass file, each line containing 5 data fields as defined
    the in the ecnf-format.txt file. Input is a list of dicts.
    """
    out = file(outfile, 'w')
    for c in curves:
        line = " ".join([c['field_label'],
                         c['N_label'],
                         c['iso_label'],
                         c['c_num'],
                         c['isomat']])
        out.write(line+"\n")
    out.close()

def rewrite_curve_file(infile, outfile, verbose=True):
    """
    Convert ideal labels for IQFs fro old-atyle N.c.d to new N.i LMFDB
    standard.  Could also be used over other fields to standardise
    labels into the canonical LMFDB ordering of ideals of the same
    norm.

    Can be used for curves* files,  curve_data* and isoclass* files.
    """
    if 'curves' in infile:
        read_file = read_curve_file
        write_file = write_curve_file
    elif 'curve_data' in infile:
        read_file = read_curvedata_file
        write_file = write_curvedata_file
    elif 'isoclass' in infile:
        read_file = read_isoclass_file
        write_file = write_isoclass_file
    else:
        print("Invalid filename {}: should be a curves or curve_data or isoclass file")
        return
    curves = read_file(infile)
    for c in curves:
        field_label = c['field_label']
        cond_label =  c['N_label']
        new_cond_label = convert_conductor_label(field_label, cond_label)
        if verbose:
            print("Conductor label {} converted to {}".format(cond_label, new_cond_label))
        c['N_label'] = new_cond_label

    write_file(curves, outfile)


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
                K = nf_lookup(field_label)
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
        using the LMFDB ordering of ideals of the same norm), then by
        isogeny class (using the key ???), then by curves in the class
        (using the key isogeny_class_key).

        By default, the letter labels for isogeny classes are
        generated on the fly, which should be fine provided that the
        input is complete (at least one curve in every isogeny class)
        for each conductor present.

        Or, the 'curves' delivered by the iterator can be tuples
        (field_label, N_label, N_def, iso_label, c_num, E) in which
        case these will be used to construct labels (as for curves
        obtained from HMFs where we want to keep the isogeny class
        label as the same as the original HMF label).

        Or, they can be generated by comparison with the associated
        newform: currently only available for Bianchi newforms over
        K=Q(sqrt(-d)) for d=1,2,3,7,11.  In these cases, since curves
        with CM rational over K are not attached to Bianchi newforms,
        the letter labels are of the form "CMa", "CMb", etc and not
        just "a", "b", etc.  """
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
                                print("processing E = %s with conductor %s, label %s-%s%s over field %s" % (list(E.ainvs()),N_def,N_label,iso_label,c_num,field_label))
                else:
                        if verbose>0:
                                print("processing E = %s..." % list(E.ainvs()))
                k = E.base_field()
                add_field(k, field_label=field_label)
                #D = Dlists[k]
                #G = Glists[k]
                used = used_curves[k]
                #isog_class_cmp = ic_cmp[k]
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
                        clist = sorted(clist0, key=isomorphism_class_key)
                        # perm[i]=j where sorted#i = unsorted#j
                        perm = dict([(i,clist0.index(Ei)) for i,Ei in enumerate(clist)])
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

                        # Q-curve? (isogeny class invariant)
                        q_curve = int(is_Q_curve(E))
                        if verbose>1:
                            if q_curve=='?':
                                print("Q-curve status not determined")
                            else:
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
				disc = cm_j_dict.get(j, 0)

                                curve_line = "%s %s %s %i %s %i %s %i %s" % (field_label, cond_label, isog, n + 1, cond_def, norm, ainv_string, disc, q_curve)
                                if verbose>1:
                                        print("curve_line: %s" % curve_line)
				tmp.append(curve_line)
			data_k[norm][hnf[1][0]][hnf[1][1]].append(tmp)

	# Sort and output the data

        #print(data)
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
				#dat.sort(cmp = isog_class_cmp)
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

# The next 2 functions are copied from lmfdb/ecnf/WebEllipticCurve.py

def ideal_from_string(K,s, IQF_format=False):
    r"""Returns the ideal of K defined by the string s.  If IQF_format is
    True, this is "[N,c,d]" with N,c,d as in a label, while otherwise
    it is of the form "[N,a,alpha]" where N is the norm, a the least
    positive integer in the ideal and alpha a second generator so that
    the ideal is (a,alpha).  alpha is a polynomial in the variable w
    which represents the generator of K (but may actially be an
    integer).  """
    #print("ideal_from_string({}) over {}".format(s,K))
    N, a, alpha = s.split(".")
    N = ZZ(N)
    a = ZZ(a)
    if IQF_format:
        d = ZZ(alpha)
        I = K.ideal(N//d, K([a, d]))
    else:
        # 'w' is used for the generator name for all fields for
        # numbers stored in the database
        alpha = alpha.encode().replace('w',str(K.gen()))
        I = K.ideal(a,K(alpha.encode()))
    if I.norm()==N:
        return I
    else:
        return "wrong" ## caller must check

def ideal_to_string(I,IQF_format=False):
    K = I.number_field()
    if IQF_format:
        a, c, d = ideal_HNF(I)
        return "[%s,%s,%s]" % (a * d, c, d)
    N = I.norm()
    a = I.smallest_integer()
    gens = I.gens_reduced()
    alpha = gens[-1]
    assert I == K.ideal(a,alpha)
    alpha = str(alpha).replace(str(K.gen()),'w')
    return ("[%s,%s,%s]" % (N,a,alpha)).replace(" ","")

def fix_conductors(infile, outfile = None, verbose=False):
    """Read curves from infile, check that their conductors as defined in
    field #5 are correct, and if not replace them.
    """
    f = outfile
    if f:
        outfile = file(f,'w')
    for L in file(infile).readlines():
        data = L.split()
        if len(data)!=13:
            print "line %s does not have 13 fields, skipping" % L
            continue
        field_label = data[0]
        K = nf_lookup(field_label)
        E = curve_from_strings(K, data[6:11])
        NE = E.conductor()
        N_def = data[4]
        N = ideal_from_string(K, N_def)
        if N!=NE:
            print("\nN_def={} gives N={}, not NE={}, replacing".format(N_def,N,NE))
            data[4] = N_def = ideal_to_string(NE).replace(" ","")
            if verbose:
                print(L[:-1])
                print(" with")
        if f:
            if verbose:
                print("writing to {}: {}".format(f,data))
            outfile.write(" ".join(data)+"\n")
        if verbose:
            print " ".join(data)
    if f:
        outfile.close()

def fix_conductor_ideals(infile, outfile = None, verbose=False):
    """Read curves from infile, replace the HNF-based format for
    conductors in field #5 with the general one.
    """
    f = outfile
    if f:
        outfile = file(f,'w')
    for L in file(infile).readlines():
        data = L.split()
        if len(data)!=13:
            print "line %s does not have 13 fields, skipping" % L
            continue
        field_label = data[0]
        K = nf_lookup(field_label)
        N_def = data[4]
        N = ideal_from_string(K, N_def, True) # input in IQF_format
        if N=="wrong":
            raise ValueError
        data[4] = N_def = ideal_to_string(N, False) # output in standard format
        if f:
            if verbose:
                print("writing to {}: {}".format(f,data))
            outfile.write(" ".join(data)+"\n")
        if verbose:
            print " ".join(data)
    if f:
        outfile.close()

def isModular(E):
    # Create a new magma instance for each curve:
    mag = Magma()
    # read in Samir Siksek's code:
    mag.eval('load "modularitycheck.m";\n')
    # Define the number field in Magma and the list of primes
    mag.eval("Qx<x> := PolynomialRing(RationalField());\n")
    K = E.base_field()
    name = K.gen()
    pol = K.defining_polynomial()
    mag.eval("Qx<x> := PolynomialRing(RationalField());\n")
    mag.eval("K<%s> := NumberField(%s);\n" % (name, pol))
    mag.eval("E := EllipticCurve(%s);\n" % list(E.ainvs()))
    mag.eval("res := isModular(E);\n")
    res = mag('res;').sage()
    mag.quit()
    return res

def check_modularity(pth,fld, verbose=False):
    f = "{}/curves.{}".format(pth,fld)
    #print("reading curves from {}".format(f))
    ntrue = 0
    nfalse = 0
    nall = 0
    for data in read_curves(f):
        field_label, N_label, N_def, iso_label, c_num, E = data
        if c_num != '1':
            continue
        c_label = N_label+"-"+iso_label+c_num
        full_label = field_label+"-"+c_label
        nall += 1
        if isModular(E):
            ntrue +=1
            if verbose: print("{} is modular".format(full_label))
        else:
            nfalse +=1
            if verbose: print("{}: could not check modularity".format(full_label))
    if ntrue==nall:
        print("All {} curves over {} were proved to be modular!".format(ntrue,field_label))
    else:
        print("Only {} out of {} curves over {} were proved to be modular!".format(ntrue,nall, field_label))

def convert_curve_file(infilename, outfilename, ncurves=0):
        r"""
        Iterator to loop through lines of a curves.* file each
        containing 13 data fields as defined the in the ecnf-format.txt file,
        writing output file in format needed by Sutherlands galdata program.
        """
        count=0
        outfile = file(outfilename, mode='w')
        for L in file(infilename).readlines():
            #sys.stdout.write(L)
            data = L.split()
            if len(data)!=13:
                print "line %s does not have 13 fields, skipping" % L
                continue
            count +=1
            if ncurves and count>ncurves:
                outfile.close()
                return

            K = nf_lookup(data[0])
            add_field(K)
            pol = pols[K]
            label = "-".join(data[:3]) + data[3]
            coeffs = ";".join(data[6:11])
            outfile.write(":".join([label,coeffs,pol])+"\n")
        outfile.close()
        return

def curve_from_data(c):
    K = nf_lookup(c['field_label'])
    return curve_from_strings(K,c['ainvs'])

def check_Q_curve_flags(filename, output=False, verbose=True):
    curves = read_curve_file(filename)
    ncurves = len(curves)
    if verbose:
        print("Read {} curves from {}".format(ncurves, filename))
    nall = 0
    nbad = 0
    nbad01 = 0
    nbad10 = 0
    cache = {} # dict of class_label:flag
    for c in curves:
        old_flag = c['q_curve_flag']
        if old_flag=='?' and not output:
            continue
        nall += 1

        class_label = "-".join([c['field_label'], c['N_label'], c['iso_label']])

        if class_label in cache:
            flag = cache[class_label]
        else:
            flag = is_Q_curve(curve_from_data(c),c['field_label'],verbose) # True, False, '?'
            cache[class_label] = flag

        if flag != '?':
            flag = str(int(flag))
        if verbose:
            print("Q-curve flag (computed): {}".format(flag))
            print("Q-curve flag (file):     {}".format(old_flag))
        #assert flag!='?'
        if old_flag!='?' and flag != old_flag:
            print("curve {}: flag on file is {} but should be {}".format(c,old_flag,flag))
            nbad += 1
        if [old_flag, flag] == ['0','1']:
            nbad01 += 1
        if [old_flag, flag] == ['1','0']:
            nbad10 += 1

        c['q_curve_flag'] = flag # a string, either '0' or '1'

    if output:
        newfilename = filename+".qc_fix"
        print("Writing revised curves file to {}".format(newfilename))
        write_curve_file(curves,newfilename)

    print("{} out of {} curves had incorrect Q-curve flag".format(nbad,nall))
    if nbad:
        print("In {} cases flag on file was 0, new flag is 1".format(nbad01))
        print("In {} cases flag on file was 1, new flag is 0".format(nbad10))
