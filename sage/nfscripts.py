# coding=utf-8
#
# Functions written to help filling in gaps in the Bianchi (imaginary
# quadratic field) curve tables, using Bianchi newforms data and Magma
# scripts.
#
from sage.all import (polygen, ZZ, QQ, Magma, magma, latex, EllipticCurve, primes, Infinity,
                      flatten, Primes, legendre_symbol, prod, RR, RealField,
                      PowerSeriesRing, O, Integer, srange, sign)
from fields import add_field, field_data, nf_lookup
from psort import nf_key, primes_of_degree_iter
from codec import ainvs_to_string, ainvs_from_string, curve_from_string, curve_from_strings, ideal_to_string, ideal_from_string, parse_point, encode_points, decode_points_one2many

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
        if I.ring()==ZZ:
            return I.gen()
        return I.norm()

def ap_list(E, Plist=None):
        r"""
        Return [a_p(E) for p in Plist].

        INPUT:

        - ``E`` - an elliptic curve defined over a number field `k`;
        - ``Plist`` - a list of primes of `k`, or None (default) in which case field_data[k]['Plist'] is used.

        OUTPUT:

        A list of a_P(E) for P in the list Plist.
        """
        if Plist is None:
            K = E.base_field()
            add_field(K)
            Plist = field_data[K]['Plist']
        return [ap(E,p) for p in Plist]

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
    G = field_data[K]['G']
    NL = G[0](N) # base-change to Galois closure
    return all([sigma(N)==NL for sigma in G.gens()])

def conj_curve(E,sigma):
    r"""
    Return the Galois conjugate elliptic curve under sigma.
    """
    return EllipticCurve([sigma(a) for a in E.ainvs()])


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
        c = [L, [L2[n] for n in sorted(L2.keys())]]
        if c:
                return c
        # For testing purposes:
        if not E1.is_isogenous(E2):
                print("Warning: curves %s and %s of conductor %s have matching L-functions\n   %s but are not isogenous!" % (E1.ainvs(),E2.ainvs(), E1.conductor(), L))
                return c

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


# Isogeny class comparison: original form, using the L-functions as
# sums over integral ideals of k.  This matches the sorting of Bianchi
# newforms.

def isog_class_cmp1(k, I, J):
    E_I = curve_from_strings(k,I[0].split()[6:11])
    E_J = curve_from_strings(k,J[0].split()[6:11])

    if not k in field_data:
        add_field(k)
    for p in field_data[k]['Plist']:
        c = int(ap(E_I, p) - ap(E_J, p))
        if c: return sign(c)

    raise NotImplementedError("Bound on primes is too small to determine...")


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

def local_data(E):
    r"""Return a local data structure, which is a list of dicts, one for each bad prime, with keys

    'p', 'normp', 'ord_cond', 'ord_disc', 'ord_den_j', 'red', 'rootno', 'kod', 'cp'

    These are all computable in Sage except for the local root number at additive primes.

    Note that The model of E might not be a global minimal model, so
    there may be one or more (in practice no more than one) entry with
    good reduction in the list. This causes no problems except that
    the bad_reduction_type is then None which cannot be converted to
    an integer.  The bad reduction types are coded as (Sage) integers
    in {-1,0,1}.

    """
    Eld = E.local_data()
    if any([ld.bad_reduction_type()==0 for ld in Eld]):
        mE = magma(E) # for local root numbers if not semistable
        mE.Conductor() # otherwise the RootNumber() function sometimes fails strangely
    def local_root_number(ldp): # ldp is a component of E.local_data()
        red_type = ldp.bad_reduction_type()
        if red_type==0: # additive reduction: call Magma
            # print("Calling Magma's RootNumber(E,P) with E = {}".format(mE))
            # print(" and P = {} = {}".format(ldp.prime(), magma(ldp.prime())))
            eps = mE.RootNumber(ldp.prime())
        elif red_type==+1:
            eps = -1
        else:  # good or non-split multiplcative reduction
            eps = +1
        return int(eps)

    E_local_data = [{'p': ideal_to_string(ld.prime()),
                   'normp': str(ld.prime().norm()),
                   'ord_cond': int(ld.conductor_valuation()),
                   'ord_disc': int(ld.discriminant_valuation()),
                   'ord_den_j': int(max(0,-(E.j_invariant().valuation(ld.prime())))),
                   'red': None if ld.bad_reduction_type() is None else int(ld.bad_reduction_type()),
                   'rootno': local_root_number(ld),
                   'kod': ld.kodaira_symbol()._pari_code(),
                   #'kod': str(latex(ld.kodaira_symbol())),

                   'cp': int(ld.tamagawa_number())}
                  for ld in Eld]
    return E_local_data, E.non_minimal_primes(), E.minimal_discriminant_ideal()

def global_period(E, scale = None, prec = None):
    r"""Return the global period of E.  This is the product over all
    infinite places v of the base field K of a local period at v,
    times a scaling factor to allow for the model not being a global
    minimal model.

    The factor at each v is E.period_lattice(v).omega().  This
    includes a factor 2 at real places where the discriminant is
    positive, i.e. the number ofc onnected components.

    The correction factor is the (rational) 12th root of the norm of
    E.discriminant()/E.minimal_discriminant_ideal().
    """
    # In Sage 9.1 there's a bug in
    # E.period_lattice(e).omega(prec=prec) for complex places where
    # the prec parameter is *not* sassed onto the period computation.
    #
    # Otherwise this would just be
    # return prod(E.period_lattice(e).omega(prec=prec) for e in K.places())

    def omega(L):
        if L.is_real():
            return L.omega(prec) if prec else L.omega()
        else:
            w1,w2 = L.basis(prec) if prec else L.basis()
            return (w1*w2.conjugate()).imag().abs()

    om = prod(omega(E.period_lattice(e)) for e in E.base_field().places())
    if scale:
        om *= scale
    return om

def extend_mwdata_one(Edata, classdata, Kfactors, magma,
                      max_sat_prime = Infinity, prec=None, verbose=False):
    r"""
    Computes analytic rank and L-value using Magma, and omega (global period).
    Computes analytic Sha (rounded).

    The prec parameter affects the precision to which the L-value and
    global period is computed.  It is bit precision.  Magma's default
    is 6dp or 20 bits for the L-value and the running time increases
    rapidly.
    """
    if prec is None:  # Magma's precision variable is decimal, 53 bits is 16 digits
        RR = RealField()
        prec = RR.precision()
        magma_prec = 16
    else:
        RR  = RealField(prec)
        # log(2)/log(10) =  0.301029995663981
        magma_prec = (prec*0.301029995663981).round()

    from fields import nf_lookup
    K = nf_lookup(Edata['field_label'])
    # We need to construct every E as a Sage EllipticCurve in
    # order to compute omega, but we only need construct it as
    # a Magma curve once per isogeny class.
    E = curve_from_string(K,Edata['ainvs'])

    # find analytic rank and L-value:

    class_label = Edata['class_label']
    if not class_label in classdata: # then we need to compute analytic rank and L-value
        mE = magma(E)
        if verbose:
            print("Calling Magma's AnalyticRank()")
        ar, lval = mE.AnalyticRank(Precision=magma_prec, nvals=2)
        if 'CM' in class_label and all(ai in QQ for ai in E.ainvs()): # avoid Magma bug
            if verbose:
                print("Special CM case: E = {}".format(E.ainvs()))
                print("AnalyticRank's ar={}, lval = {}".format(ar,lval))
            ar *= 2
            old_lval = lval
            lval = mE.LSeries().Evaluate(1, Derivative=ar) / magma.Factorial(ar)
            if verbose:
                print("ar doubled to {}, lval recomputed to {}".format(ar,lval))
                print(" (compare square of old lval:       {})".format(old_lval**2))
        lval = RR(lval)
        ar = int(ar)
        classdata[class_label] = (ar,lval)
    else:
        ar, lval = classdata[class_label]
    Edata['analytic_rank'] = ar
    Edata['Lvalue'] = lval
    if verbose:
        print("analytic rank = {}\nL-value = {}".format(ar,lval))

    # recompute regulator.  Original heights were computed
    # before fixing Sage's height function precision issues
    # properly.

    gens = [E(parse_point(K,P)) for P in Edata['gens']]
    ngens = len(gens)
    if verbose:
        print("gens = {}".format(gens))

    if max_sat_prime and ngens:
        if max_sat_prime==Infinity:
            try:
                new_gens, index, new_reg = E.saturation(gens, verbose=verbose)
            except ValueError:
                print("Warning: unable to compute saturation index bound, using 100")
                new_gens, index, new_reg = E.saturation(gens, max_prime=100, verbose=verbose)
        else:
            new_gens, index, new_reg = E.saturation(gens, max_prime=max_sat_prime, verbose=verbose)
        if index>1:
            print("Original gens were not saturated, index = {} (using max_prime {})".format(index,max_sat_prime))
            gens = new_gens
            Edata['gens'] = decode_points_one2many(encode_points(gens)) # list of strings
        else:
            if verbose:
                print("gens are saturated at primes up to {}".format(max_sat_prime))

    heights = [P.height(precision=prec) for P in gens]
    Edata['heights'] = str(heights).replace(" ","")
    if verbose:
        print("heights = {}".format(heights))
    reg = E.regulator_of_points(gens, precision=prec)
    Edata['reg'] = str(reg) if ar else '1'
    if verbose:
        print("regulator (of known points) = {}".format(reg))

    # allow for the scaling in the Neron-Tate height
    # pairing: for BSD we need non-normalised heights and
    # normalization divides every height by K.degree(), so
    # the regulator we need has to be multiplied by
    # K.degree()**rank.
    if len(gens) == ar:
        NTreg = reg * K.absolute_degree()**ar
        if verbose:
            print("Neron-Tate regulator = {}".format(NTreg))
    else:
        NTreg = None

    # compute omega

    # find scaling factor in case we don't have a global minimal model
    minDnorm = ZZ(Edata['minD'][1:].split(",")[0]).abs()
    modelDnorm = E.discriminant().norm().abs()
    fac = (modelDnorm/minDnorm).nth_root(12) # will be exact
    if fac!=1 and verbose:
        print("Not a global minimal model")
        print("Scaling factor = {}".format(fac))
    Edata['omega'] = omega = global_period(E, fac, prec=prec)
    if verbose:
        print("omega = {}".format(omega))

    T = E.torsion_subgroup()
    nt = T.order()
    if nt != Edata['torsion_order']:
        print("{}: torsion order is {}, not {} as on file; updating data".format(Edata['label'],nt,Edata['torsion_order']))
        Edata['torsion_order'] = nt
        Edata['torsion_structure'] = list(T.invariants())
        tgens = [P.element() for P in T.gens()]
        Edata['torsion_gens'] = decode_points_one2many(encode_points(tgens)) # list of strings
    if verbose:
        print("Torsion order = {} (checked)".format(nt))

    tamagawa_product = Edata['tamagawa_product']
    if verbose:
        print("Tamagawa product = {}".format(tamagawa_product))

    if NTreg:
        Rsha = lval * nt**2  / (NTreg * tamagawa_product * omega)

        if not K in Kfactors:
            Kfactors[K] = RR(K.discriminant().abs()).sqrt() / 2**(K.signature()[1])
        if verbose:
            print("Field factor = {}".format(Kfactors[K]))

        Rsha *= Kfactors[K]
        Edata['sha'] = sha = Rsha.round()
        if verbose:
            print("Approximate analytic Sha = {}, rounds to {}".format(Rsha, sha))
        if sha==0 or (sha-Rsha).abs()>0.0001 or not ZZ(sha).is_square():
            if not verbose:
                print("Approximate analytic Sha = {}, rounds to {}".format(Rsha, sha))
            print("****************************Not good! 0 or non-square or not close to a positive integer!")

    else:
        if verbose:
            print("Unable to compute regulator or analytic Sha, since analytic rank = {} but we only have {} generators".format(ar, Edata['ngens']))
        Edata['sha'] = None
    return Edata

def latex_equation(ainvs):
    a1, a2, a3, a4, a6 = ainvs

    def co(coeff):
        pol = coeff.polynomial()
        mons = pol.monomials()
        n = len(mons)
        if n==0:
            return ""
        if n>1:
            return r"+\left({}\right)".format(latex(coeff))
        # now we have a numerical coefficient times a power of the generator
        if coeff == 1:
            return "+"
        if coeff == -1:
            return "-"
        s = "+" if pol.monomial_coefficient(mons[0]) > 0 else ""
        return "{}{}".format(s, latex(coeff))

    def term(coeff, mon):
        if not coeff:
            return ""
        if not mon:
            return "+{}".format(latex(coeff)).replace("+-","-")
        return "{}{}".format(co(coeff), mon)

    return ''.join([r'y^2',
                    term(a1,'xy'),
                    term(a3,'y'),
                    '=x^3',
                    term(a2,'x^2'),
                    term(a4,'x'),
                    term(a6,''),
                    r'']).replace(" ","")

def simplify_one_ideal_string(K, Istring):
    """Given an ideal string in the form "[N,a,alpha]" representing an
    ideal of K with norm N, return a string of the form "(g)" or "(g1,
    g2)" defining the same ideal with the minimal number of
    generators.
    """
    return str(ideal_from_string(K, Istring).gens_reduced()).replace(",)",")").replace(" ","")

def ec_disc(ainvs):
    """
    Return disciminant of a Weierstrass equation from its list of a-invariants.
    Avoids constructing the EllipticCurve.
    """
    a1, a2, a3, a4, a6 = ainvs
    b2 = a1*a1 + 4*a2
    b4 = a3*a1 + 2*a4
    b6 = a3*a3 + 4*a6
    c4 = b2*b2 - 24*b4
    c6 = -b2*b2*b2 + 36*b2*b4 - 216*b6
    return (c4*c4*c4 - c6*c6) / 1728

def reduce_mod_units(a):
    """
    Return u*a for a unit u such that u*a is reduced.
    """
    K = a.parent()
    if a.norm().abs()==1:
        return K(1)
    r1, r2 = K.signature()
    if r1 + r2 == 1:  # unit rank is 0
        return a

    from sage.modules.all import vector
    from sage.matrix.all import Matrix

    prec = 1000  # lower precision works badly!

    embs = K.places(prec=prec)
    aconjs = [e(a) for e in embs]
    degs = [1]*r1 + [2]*r2
    v = vector([aa.abs().log()*d for aa,d in zip(aconjs,degs)])

    fu = K.units()
    U = Matrix([[e(u).abs().log()*d for d,e in zip(degs,embs)] for u in fu])
    A = U*U.transpose()
    Ainv = A.inverse()

    exponents = [e.round() for e in -Ainv*U*v]
    u = prod([uj**ej for uj,ej in zip(fu,exponents)])
    return a*u

def simplify_ideal_strings(K, record):
    """Convert ideal strings from long form [N,a,alpha] to 1-generator
    (gen) or 2-generators (gen1,gen2).
    """
    if not K:
        K = nf_lookup(record['field_label'])
    record['conductor_ideal'] = simplify_one_ideal_string(K, record['conductor_ideal'])
    record['minD'] = simplify_one_ideal_string(K, record['minD'])
    #print(record['bad_primes'])
    record['bad_primes'] = [simplify_one_ideal_string(K, p) for p in record['bad_primes']]
    #print(record['bad_primes'])
    #print(record['non_min_p'])
    record['non_min_p'] = [simplify_one_ideal_string(K, p) for p in record['non_min_p']]
    for ld in record['local_data']:
        ld['p'] = simplify_one_ideal_string(K, ld['p'])
    # avoid constructing the curve! D is the principal ideal generated
    # by the discriminant.
    ainvs = ainvs_from_string(K, record['ainvs'])
    D = ec_disc(ainvs)
    Dnorm = D.norm()
    if Dnorm==1:
        D = 1
    else:
        Dred = reduce_mod_units(D)
        if len(str(Dred)) < len(str(D)):
           D = Dred
    record['disc'] = '({})'.format(D).replace(" ","")
    record['normdisc'] = Dnorm
    return record

# The following was used (with fix_models_field() in files.py) to
# adjust the data for any field where the stored models were not
# minimal.

def fix_model(K, record, verbose=False):
    """
    Change to a (semi-)global minimal model.  Keys affected are:

    ainvs, equation (in curves.*)
    local_data, non_min_p (in local_data.*)
    (in local_data there may be fewer primes, but the data at each prime should be unchanged)
    gens, torsion_gens (in mwdata.*)
    """
    if verbose:
        print("fixing data for {}".format(record['label']))
    ainvs = ainvs_from_string(K, record['ainvs'])
    E = EllipticCurve(ainvs)
    Emin = E.global_minimal_model(semi_global=True)
    iso = E.isomorphism_to(Emin)

    # change ainvs and equation:
    if verbose:
        print("old ainvs: {}".format(record['ainvs']))
    record['ainvs'] = ainvs_to_string(Emin.ainvs())
    if verbose:
        print("new ainvs: {}".format(record['ainvs']))
        print("old equation: {}".format(record['equation']))
    record['equation'] = latex_equation(Emin.ainvs())
    if verbose:
        print("new equation: {}".format(record['equation']))

    # map gens and torsion gens:
    if verbose:
        print("old gens: {}".format(record['gens']))
    new_gens = [iso(E(parse_point(K, P))) for P in record['gens']]
    record['gens'] = decode_points_one2many(encode_points(new_gens))
    if verbose:
        print("new gens: {}".format(record['gens']))
        print("old torsion_gens: {}".format(record['torsion_gens']))
    new_torsion_gens = [iso(E(parse_point(K, P))) for P in record['torsion_gens']]
    record['torsion_gens'] = decode_points_one2many(encode_points(new_torsion_gens))
    if verbose:
        print("new torsion_gens: {}".format(record['torsion_gens']))

    # new local_data:
    ld, non_min_p, minD = local_data(Emin)

    if verbose:
        print("old local_data: {}".format(record['local_data']))
        print("new local_data: {}".format(ld))
        print("old non_min_p: {}".format(record['non_min_p']))
        print("new non_min_p: {}".format(non_min_p))
        print("old minD: {}".format(record['minD']))
        print("new minD: {}".format(minD))
    record['local_data'] = ld
    record['non_min_p'] = [ideal_to_string(P) for P in non_min_p]
    record['minD'] = ideal_to_string(minD)

    #record = simplify_ideal_strings(K, record)
    return record
