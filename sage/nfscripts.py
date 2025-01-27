# coding=utf-8
#
# Functions for computing elliptic curve data
#
from os import path
from sage.all import (polygen, ZZ, QQ, latex,
                      EllipticCurve, primes, flatten, Primes,
                      legendre_symbol, prod, RealField,
                      PowerSeriesRing, O, Integer, srange, sign, copy)

from files import ECNF_DIR

from fields import (add_field, field_data, cm_j_dict, get_field_label, get_field_type_from_label, nf_lookup, subfield_labels)

from codec import (ainvs_to_string, ainvs_from_string,
                   curve_from_string, curve_from_strings,
                   ideal_to_string, ideal_from_string, parse_point,
                   parse_curves_line, encode_points, decode_points_one2many)

from magma import get_magma

def get_Q_label(EQ):
    try:
        return EQ.label()
    except LookupError:
        return "{}.?.?".format(EQ.conductor())

# Given an elliptic curve, find it in the (text file) database and
# return its label.  Note that we can filter by field_label and
# conductor_norm, but not by conductor_label, since the curves in the
# database over totally real fields get their conductor label from the
# associated HMF label, which comes from Magma and is *not* in general
# the standard LMFDB label foe the ideal.
def get_label(E):
    K = E.base_field()
    if K.degree()==1:
        return get_Q_label(E)
    field_label = get_field_label(K)
    field_type = get_field_type_from_label(field_label)
    cond = E.conductor()
    cond_norm = cond.norm()
    curves_file = path.join(ECNF_DIR, field_type, "curves."+field_label)
    #print("Looking for {}".format(E.ainvs()))
    #print("curves file: {}".format(curves_file))
    with open(curves_file) as infile:
        for line in infile.readlines():
            label, data = parse_curves_line(line)
            if data['conductor_norm'] != cond_norm:
                continue
            #print(" - testing against {} = {}".format(label, data['ainvs']))
            E1 = curve_from_string(K, data['ainvs'])
            if E.is_isomorphic(E1):
                #print("Success! Label is {}".format(label))
                return label
    from psort import ideal_label
    cond_label = ideal_label(cond)
    lab = "-".join([field_label, cond_label, "?", "?"])
    #print("*************** {} not found, returning {}".format(E.ainvs(), lab))
    return lab

# return a list of labels of elliptic curves over proper subfields
# whose base-change to this curve's base field is this curve.  When
# one or more of the curves over subfields is not in the database, the
# label will have '?' for the isogeny class label and the number
# within the class (though the latter could be computed).
def get_base_change_labels(E):
    K = E.base_field()
    degK = K.degree()
    fromQ = [get_Q_label(E0) for E0 in E.descend_to(QQ)]
    # if fromQ:
    #     return fromQ
    if degK==1 or degK.is_prime():
        return fromQ
    field_label = get_field_label(K, exact=False) # allow isomorphism
    sub_labels = subfield_labels(field_label)
    return fromQ+flatten([[get_label(E0) for E0 in E.descend_to(nf_lookup(lab))] for lab in sub_labels])

def get_all_base_change(outfile=None, ncurves=0):
    from fields import nf_table
    if (outfile):
        out = open(outfile, 'w')
        out.write("label|base_change\n")
        out.write("text|jsonb\n\n")
    bc_dict = {}
    for field_label, K in nf_table.items():
        if K.degree().is_prime():
            continue # silently skip fields of degree 2,3,5
        sf = subfield_labels(field_label)
        if not sf:
            print("Skipping {} as it has no nontrivial proper subfields".format(field_label))
            continue
        with open(path.join(ECNF_DIR, get_field_type_from_label(field_label), "curves.{}".format(field_label))) as infile:
            print("Processing {} with nontrivial proper subfields {}".format(field_label, sf))
            nc = 0
            for L in infile.readlines():
                label, data = parse_curves_line(L)
                E = curve_from_string(K, data['ainvs'])
                bc_labels = get_base_change_labels(E)
                if bc_labels:
                  bc_dict[label] = bc_labels
                  print("Curve {} is bc of {}".format(label, bc_labels))
                  if (outfile):
                      out.write("{}|{}\n".format(label, bc_labels))
                nc +=1
                if nc%1000==0:
                  print("checked {} curves over {} so far".format(nc, field_label))
                if ncurves and nc==ncurves:
                    if (outfile):
                        out.close()
                    return bc_dict

            print("Finished checking {} curves over {}".format(nc, field_label))
    if (outfile):
        out.close()
    return bc_dict

def torsion_data(E):
    T = E.torsion_subgroup()
    tdata = {}
    tdata['torsion_structure'] = ts = list(T.invariants())
    tdata['torsion_gens'] = tg = [P.element() for P in T.smith_form_gens()]
    tdata['torsion_order'] = T.order()
    if len(ts) == 2:
        assert ts[1]%ts[0] == 0
        assert tg[0].order() == ts[0]
        assert tg[1].order() == ts[1]
    return tdata

def ap(E, p):
    r""" Return a_p(E).

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
    if I.ring() == ZZ:
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
    return [ap(E, p) for p in Plist]

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
    if K is QQ:
        return True
    add_field(K, field_label=field_label)
    G = field_data[K]['G']
    NL = G[0](N) # base-change to Galois closure
    return all([sigma(N) == NL for sigma in G.gens()])

def conj_curve(E, sigma):
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
            yield p

def cmj_key(E):
    r""" Key to compare curves with non-rational CM which are quadratic
    twists over the CM field.  This will be called on lots of curves
    for which this tie-break comparison is not needed, so we return 0
    instantly when we know that is the case.
    """
    from psort import primes_of_degree_iter
    if (not E.has_cm()) or E.has_rational_cm():
        return 0
    d = E.cm_discriminant()
    K = E.base_field()
    deg = K.absolute_degree()
    D = 1 if deg == 1 else ZZ(K.defining_polynomial().discriminant())
    j = E.j_invariant()
    c4, c6 = E.c_invariants()
    jj, c, w = (j, c4, 4) if j == 1728 else (j-1728, c6, 6)
    NN = E.conductor() if deg == 1 else E.conductor().norm()
    bad = 6*d*D*NN

    # Get the first degree 1 prime P, dividing a prime p not dividng
    # bad for which d is a quadratic non-residue, such that j-1728 (or
    # j when j=1728) is a P-unit:
    ptest = lambda p: not p.divides(bad) and legendre_symbol(d, p) == -1
    if deg == 1:
        it = primes_iter_Q(ptest)
    else:
        it = primes_of_degree_iter(K, deg=1, condition=ptest)
    P = next(it)
    while jj.valuation(P) != 0:
        P = next(it)
    p = P if deg == 1 else P.smallest_integer() # = residue characteristic
    print("E = {} with j = {}: tie-break prime P = {} above p = {}".format(E.ainvs(), j, P, p))

    # The key is now (c6|p) (or (c4|p) if c6=0) with c4, c6 from the
    # P-minimal model.  Although E has good reduction at P the model
    # may not be minimal, and some adjustment is necessary:
    k = c.valuation(P)
    if k > 0:
        w = ZZ(w)
        assert w.divides(k)
        pi = K.uniformizer(P, others='negative')
        c = c/pi**(k//w) # still integral everywhere
    return legendre_symbol(K.residue_field(P)(c), p)

def isomorphism_class_key_j(E):
    """For isogenous curves, first sort by CM-discriminant, then by
    j-invariant, then (only necessary when E has potential CM) the
    tie-break.

    """
    from psort import nf_key
    return (int(E.has_rational_cm() and -E.cm_discriminant()),
            nf_key(E.j_invariant()),
            cmj_key(E))

isomorphism_class_key = isomorphism_class_key_j

def Euler_polynomial(E, P):
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
        return 1 - EP.bad_reduction_type()*polygen(ZZ)

def rational_Euler_polynomial(E, p):
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
    return prod([Euler_polynomial(E, P)(x^P.residue_class_degree())
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
    from sage.all import RR
    # maxexp(p) = max{i: p^i <= nmax}
    lognmax = RR(nmax).log()
    maxexp = lambda p: (lognmax/RR(p).log()).floor()

    polydata = [(p, maxexp(p), rational_Euler_polynomial(E, p))
                for p in primes(nmax+1)]
    t = PowerSeriesRing(ZZ, 't').gen()
    c = {}
    for p, e, pol in polydata:
        s = (1/pol(t) + O(t**nmax)).dict()
        cp = dict([(p^i, s.get(i, 0)) for i in range(1, e+1)])
        c.update(cp)

    # so far, c[n] is defined for n=p^i with 1<n<=nmax, but only when c[n]!=0
    c[1] = Integer(1)
    if prime_powers_only:
        return c

    for n in srange(2, nmax+1):
        if n not in c: # we do not yet know c[n]
            nf = n.factor()
            assert len(nf) > 1
            n1 = nf[0][0]**nf[0][1] # p**e
            n2 = n//n1 # n2<n so we have c[n2]
            c[n] = c[n1]*c[n2]
    return c

def curve_cmp_via_L(E1, E2, nmax=100):
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
    L1 = rational_L_coefficients(E1, nmax)
    L2 = rational_L_coefficients(E2, nmax)
    L = [L1[n] for n in sorted(L1.keys())]
    c = [L, [L2[n] for n in sorted(L2.keys())]]
    if c:
        return c
    # For testing purposes:
    if not E1.is_isogenous(E2):
        print("Warning: curves %s and %s of conductor %s have matching L-functions\n   %s but are not isogenous!" % (E1.ainvs(), E2.ainvs(), E1.conductor(), L))
        return c

# Isogeny class comparison: experimental for, based on comparison
# between the L-functions as rational Dirichlet series (indexed by
# postive integers, hence no choices needed!)  implemented at ICTP
# September 2014 by Angelos Koutsianas, Alejandro Argaez, Daniel
# Kohen, and revised here by John Cremona.  NB This has been
# *discarded* since Galois conjugate classes have the same rational
# L-function and would need a tie-break anyway.

def isog_class_cmp2(k, I, J):
    E1 = curve_from_strings(k, I[0].split()[6:11])
    E2 = curve_from_strings(k, J[0].split()[6:11])
    return curve_cmp_via_L(E1, E2)

# Isogeny class comparison: original form, using the L-functions as
# sums over integral ideals of k.  This matches the sorting of Bianchi
# newforms.

def isog_class_cmp1(k, I, J):
    E_I = curve_from_strings(k, I[0].split()[6:11])
    E_J = curve_from_strings(k, J[0].split()[6:11])

    if k not in field_data:
        add_field(k)
    for p in field_data[k]['Plist']:
        c = int(ap(E_I, p) - ap(E_J, p))
        if c:
            return sign(c)

    raise NotImplementedError("Bound on primes is too small to determine...")


def isModular_magma(E):
    # Create a new magma instance for each curve:
    mag = get_magma()
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
    return res

def local_data(E):
    r"""Return disct containing local data of E:

    'local_data': a list of dicts, one for each bad prime, with keys
    'p', 'normp', 'ord_cond', 'ord_disc', 'ord_den_j', 'red',
    'rootno', 'kod', 'cp'.  All are integers except 'p' which is a
    simplified ideal string defining the prime.

    'non_min_p': a list of non-minimal primes (as simplified ideal strings)

    'minD':  the minimal discriminant ideal (as a simplified ideal string)

    'bad_primes': a list of the primes of bad reduction (as simplified
    ideal strings). Note that this may be shorter than 'local_data' if
    there are primes of good reduction at which the model is
    non-minimal.

    These are all computable in Sage except for the local root number at additive primes.

    Note that the model of E might not be a global minimal model, so
    there may be one or more (in practice no more than one) entry with
    good reduction in the list. This causes no problems except that
    the bad_reduction_type is then None which cannot be converted to
    an integer.  The bad reduction types are coded as (Sage) integers
    in {-1,0,1}.

    """
    magma = get_magma()
    Eld = E.local_data()
    if any([ld.bad_reduction_type() == 0 for ld in Eld]):
        mE = magma(E) # for local root numbers if not semistable
        mE.Conductor() # otherwise the RootNumber() function sometimes fails strangely
    def local_root_number(ldp): # ldp is a component of E.local_data()
        red_type = ldp.bad_reduction_type()
        if red_type == 0: # additive reduction: call Magma
            # print("Calling Magma's RootNumber(E,P) with E = {}".format(mE))
            # print(" and P = {} = {}".format(ldp.prime(), magma(ldp.prime())))
            eps = mE.RootNumber(ldp.prime())
        elif red_type == +1:
            eps = -1
        else:  # good or non-split multiplcative reduction
            eps = +1
        return int(eps)

    E_local_data = [{'p': ideal_to_string(ld.prime()),
                     'normp': str(ld.prime().norm()),
                     'ord_cond': int(ld.conductor_valuation()),
                     'ord_disc': int(ld.discriminant_valuation()),
                     'ord_den_j': int(max(0, -(E.j_invariant().valuation(ld.prime())))),
                     'red': None if ld.bad_reduction_type() is None else int(ld.bad_reduction_type()),
                     'rootno': local_root_number(ld),
                     'kod': ld.kodaira_symbol()._pari_code(),
                     'cp': int(ld.tamagawa_number())}
                    for ld in Eld]
    return {'local_data': E_local_data,
            'non_min_p': [ideal_to_string(P) for P in E.non_minimal_primes()],
            'minD': ideal_to_string(E.minimal_discriminant_ideal()),
            'bad_primes': [ldp['p'] for ldp in E_local_data if ldp['ord_cond']],
           }

def global_period(E, scale=None, prec=None):
    r"""Return the global period of E.  This is the product over all
    infinite places v of the base field K of a local period at v,
    times a scaling factor to allow for the model not being a global
    minimal model.

    The factor at each v is E.period_lattice(v).omega().  This
    includes a factor 2 at real places where the discriminant is
    positive, i.e. the number of connected components.

    The correction factor is the (rational) 12th root of the norm of
    E.discriminant()/E.minimal_discriminant_ideal().
    """

    # E.period_lattice(e).omega() now has a parameter bsd_normalise, False by default
    # and ignored for real places, which multiplies by 2 for complex places.

    def omega(L):
        return L.omega(prec, bsd_normalise=True)

    om = prod(omega(E.period_lattice(e)) for e in E.base_field().places())
    if scale:
        om *= scale
    return om

def analytic_rank_and_lvalue_magma(E, prec=None, verbose=False):
    from magma import get_magma
    mag = get_magma()
    if prec is None:  # Magma's precision variable is decimal, 53 bits is 16 digits
        R = RealField()
        magma_prec = 6
    else:
        R = RealField(prec)
        magma_prec = R(prec*0.301029995663981).round()
    if verbose:
        print(f"Calling Magma's AnalyticRank() with decimal precision {magma_prec}")
    ar, lval = mag(E).AnalyticRank(Precision=magma_prec, nvals=2)
    if verbose:
        print(f"Magma returns {ar=}, {lval=}")
    return int(ar), R(lval)

def analytic_rank_and_lvalue_pari(E, prec=None, rmax=4, verbose=False):
    from sage.all import pari
    if prec is None:  # Magma's precision variable is decimal, 53 bits is 16 digits
        R = RealField()
        prec = R.precision()
    else:
        R = RealField(prec)
    eps = R(2)**(1-prec)
    rmax1 = rmax+1
    pE = pari(E)
    keep_prec = pari.default('realbitprecision')
    pari.set_real_precision_bits(prec)
    assert pari.default('realbitprecision')==prec
    if verbose:
        print(f"Calling Pari's lfun at 1+x+O(x^{rmax1}), with bit precision {prec}")
    v = pE.lfun(f'1+x+O(x^{rmax1})', precision=prec)
    vc = [R(v.polcoef(i)) for i in range(rmax+1)]
    ar = next(i for i in range(rmax+1) if vc[i].abs() > eps)
    lval = vc[ar]
    if verbose:
        print(f"Pari returns coefficients {vc}\n so {ar=}, {lval=}")
    pari.set_real_precision_bits(keep_prec)
    return ar, lval

def analytic_rank_and_lvalue(E, prec=128, backend='Magma', verbose=False):
    if backend=='Magma':
        return analytic_rank_and_lvalue_magma(E, prec=prec, verbose=verbose)
    if backend=='pari':
        return analytic_rank_and_lvalue_pari(E, prec=prec, verbose=verbose)
    raise ValueError("backend value must be 'Magma' or 'pari'")


def extend_mwdata_one(Edata, classdata, Kfactors,
                      max_sat_prime=None, prec=None, backend='Magma', verbose=False):
    r"""Computes analytic rank and L-value using Magma or pari, and omega (global
    period), heights and regulator using Sage.  Computes analytic Sha
    (rounded).

    The prec parameter controls the precision to which the heights,
    regulator and global period is computed.  It is bit precision.
    """
    from fields import nf_lookup
    K = nf_lookup(Edata['field_label'])

    if prec is None:
        R = RealField()
        prec = R.precision()
    else:
        R = RealField(prec)

    # We need to construct every E as a Sage EllipticCurve in
    # order to compute omega, but we only need construct it as
    # a Magma/Pari curve once per isogeny class for the a.r. and Lvalue.
    E = curve_from_string(K, Edata['ainvs'])

    # find analytic rank and L-value (if degree<6):

    class_label = Edata['class_label']
    if class_label not in classdata: # then we need to compute analytic rank and L-value
        ar, lval = analytic_rank_and_lvalue(E, prec=prec, backend=backend, verbose=verbose)
        classdata[class_label] = (ar, lval)
    else:
        if verbose:
            print("ar and Lval already computed for isogeny class {}: {}".format(class_label, classdata[class_label]))
        ar, lval = classdata[class_label]
    Edata['analytic_rank'] = ar
    Edata['Lvalue'] = lval
    if verbose:
        print("analytic rank = {}\nL-value = {}".format(ar, lval))

    # recompute regulator.  Original heights were computed
    # before fixing Sage's height function precision issues
    # properly.

    gens = [E(parse_point(K, P)) for P in Edata['gens']]
    ngens = len(gens)
    if verbose:
        print("gens = {}".format(gens))

    if ngens:
        if max_sat_prime:
            new_gens, index, _ = E.saturation(gens, max_prime=max_sat_prime, verbose=0)
        else:
            new_gens, index, _ = E.saturation(gens, verbose=0)
        if index > 1:
            print(f"Original gens were not saturated, {index = }")
            if max_sat_prime:
                print(f" (using {max_sat_prime = })")
            gens = new_gens
            Edata['gens'] = decode_points_one2many(encode_points(gens)) # list of strings
        else:
            if verbose:
                if max_sat_prime:
                    print(f"gens are saturated at primes up to {max_sat_prime}")
                else:
                    print("gens are saturated at all primes")

    heights = [P.height(precision=prec) for P in gens]
    Edata['heights'] = str(heights).replace(" ", "")
    if verbose:
        print("heights = {}".format(heights))
    reg = E.regulator_of_points(gens, precision=prec) if gens else 1
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
    from fields import special_names
    minD = Edata['minD'].replace('w', special_names.get(Edata['field_label'], 'a'))
    minD = K.ideal([K(s) for s in minD[1:-1].split(",")])
    minDnorm = minD.norm()
    modelDnorm = E.discriminant().norm().abs()
    fac = (modelDnorm/minDnorm).nth_root(12) # will be exact
    if fac != 1 and verbose:
        print("Not a global minimal model")
        print("Scaling factor = {}".format(fac))
    Edata['omega'] = omega = global_period(E, fac, prec=prec)
    if verbose:
        print("omega = {}".format(omega))

    tdata = torsion_data(E)
    nt = tdata['torsion_order']
    if nt != Edata['torsion_order']:
        print("{}: torsion order is {}, not {} as on file; updating data".format(Edata['label'], nt, Edata['torsion_order']))
        Edata['torsion_order'] = nt
        Edata['torsion_structure'] = tdata['torsion_structure']
        tgens = tdata['torsion_gens']
        Edata['torsion_gens'] = decode_points_one2many(encode_points(tgens)) # list of strings
    if verbose:
        print("Torsion order = {} (checked)".format(nt))

    tamagawa_product = Edata['tamagawa_product']
    if verbose:
        print("Tamagawa product = {}".format(tamagawa_product))

    if NTreg:
        if K not in Kfactors:
            Kfactors[K] = R(K.discriminant().abs()).sqrt()
        if verbose:
            print("Field factor = {}".format(Kfactors[K]))

        # NB we may have computed lval to lower precision than NTreg and omega

        Rsha = Kfactors[K] * lval * nt**2  / R(NTreg * tamagawa_product * omega)
        Edata['sha'] = sha = Rsha.round()
        if verbose:
            print("Approximate analytic Sha = {}, rounds to {}".format(Rsha, sha))
        if sha == 0 or (sha-Rsha).abs() > 0.0001 or not ZZ(sha).is_square():
            if not verbose:
                print("Approximate analytic Sha = {}, rounds to {}".format(Rsha, sha))
            print("****************************Not good! 0 or non-square or not close to a positive integer!")
            if not ZZ(sha).is_square():
                print("Not keeping this Sha value")
                sha = None

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
        if n == 0:
            return ""
        if n > 1:
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
            return "+{}".format(latex(coeff)).replace("+-", "-")
        return "{}{}".format(co(coeff), mon)

    return ''.join([r'y^2',
                    term(a1, 'xy'),
                    term(a3, 'y'),
                    '=x^3',
                    term(a2, 'x^2'),
                    term(a4, 'x'),
                    term(a6, ''),
                    r'']).replace(" ", "")

def simplify_one_ideal_string(K, Istring):
    """Given an ideal string in the form "[N,a,alpha]" representing an
    ideal of K with norm N, return a string of the form "(g)" or "(g1,
    g2)" defining the same ideal with the minimal number of
    generators.
    """
    return ideal_to_string(ideal_from_string(K, Istring))

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
    if a.norm().abs() == 1:
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
    v = vector([aa.abs().log()*d for aa, d in zip(aconjs, degs)])

    fu = K.units()
    U = Matrix([[e(u).abs().log()*d for d, e in zip(degs, embs)] for u in fu])
    A = U*U.transpose()
    Ainv = A.inverse()

    exponents = [e.round() for e in -Ainv*U*v]
    u = prod([uj**ej for uj, ej in zip(fu, exponents)])
    return a*u

def simplified_gens(I):
    """
    return simplified generators for the ideal I
    """
    return [reduce_mod_units(a) for a in I.gens_reduced()]

def simplify_ideal_strings(K, record):
    """Convert ideal strings from long form [N,a,alpha] to 1-generator
    (gen) or 2-generators (gen1,gen2).
    """
    from fields import nf_lookup
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
    if Dnorm == 1:
        D = 1
    else:
        Dred = reduce_mod_units(D)
        if len(str(Dred)) < len(str(D)):
            D = Dred
    record['disc'] = '({})'.format(D).replace(" ", "")
    record['normdisc'] = Dnorm
    return record

def make_isogeny_class(curve, verbose=False, prec=None, backend='Magma'):
    """curve is a dict with keys

    field_label, conductor_norm, conductor_label, conductor_ideal, iso_label, ainvs

    defining a single curve.

    Returns a dict with keys curve labels, one for each curve in the
    isogeny class, and values complete curve records.

    """
    from mwinfo import compute_mwdata
    from trace_hash import TraceHash
    from Qcurves import is_Q_curve
    from galrep import get_galrep_data

    class_label = "-".join([curve[k] for k in ['field_label', 'conductor_label', 'iso_label']])
    if True:
        print("Processing isogeny class {}".format(class_label))

    ainvs = curve['ainvs']
    K = ainvs[0].parent()
    add_field(K)
    E = EllipticCurve(K, ainvs)

    from sage.schemes.elliptic_curves.isogeny_class import IsogenyClass_EC_NumberField
    Cl = IsogenyClass_EC_NumberField(E, reducible_primes=None, algorithm='Billerey', minimal_models=True)
    #Cl = E.isogeny_class()
    clist0 = [minimal_model(C) for C in Cl.curves]
    mat0 = Cl.matrix()
    clist = sorted(clist0, key=isomorphism_class_key)
    nc = len(clist)
    # perm[i]=j where sorted#i = unsorted#j
    perm = dict([(i, clist0.index(Ei)) for i, Ei in enumerate(clist)])
    mat = [[mat0[perm[i], perm[j]] for j in range(nc)] for i in range(nc)]

    curve['isogeny_matrix'] = mat
    curve['q_curve'] = is_Q_curve(E) # isogeny-invariant
    curve['trace_hash'] = TraceHash(E)

    # this adds 'mwdata' to curve, being a list of dicts with keys
    #'rank_bounds', 'rank', 'gens', 'hts', 'reg', 'ntors',
    #'torstruct', 'tgens'.  Some of this is isogeny-invariant, the
    #rest will get copied across to the individual curve records

    curve['curves'] = clist
    mwdata = compute_mwdata(curve, backend=backend, verbose=verbose, prec=prec)
    if verbose > 1:
        print("Computed mwdata")
    # copy class-invariant data
    curve['analytic_rank'] = mwdata[0].get('analytic_rank', None)
    curve['rank_bounds'] = mwdata[0].get('rank_bounds', None)
    curve['rank'] = mwdata[0].get('rank', None)

    # fill in data for individual curves

    data = {} # will hold all data to be returned
    for i, c in enumerate(clist):
        record = copy(curve)

        record['number'] = i + 1 # count from 1 here
        record['label'] = label = class_label + str(record['number'])
        record['ainvs'] = ainvs = c.ainvs()
        record['jinv'] = j = c.j_invariant()
        D = c.discriminant()
        Dnorm = D.norm()
        if Dnorm == 1:
            D = 1
        else:
            if D in ZZ: # fairly common special case
                D = ZZ(D).abs()
            else:
                Dred = reduce_mod_units(D)
                if len(str(Dred)) < len(str(D)):
                    D = Dred
        record['disc'] = '({})'.format(D).replace(" ", "")
        record['normdisc'] = Dnorm
        record['equation'] = latex_equation(ainvs)
        record['cm'] = cm_j_dict.get(j, 0)
        if verbose:
            print("j-invariant = {}, CM = {}".format(j, record['cm']))

        # The following will fail for curves which are
        # base-changes of curves over Q with large conductor (and
        # hence no Cremona label)
        EEQ = c.descend_to(QQ)
        if EEQ:
            if verbose:
                print("curves     over Q: {}".format([E.ainvs() for E in EEQ]))
                print("conductors over Q: {}".format([E.conductor() for E in EEQ]))
        record['base_change'] = [get_Q_label(cQ) for cQ in c.descend_to(QQ)]

        # local_data (keys 'local_data', 'non_min_p', 'minD', 'bad_primes'):
        ld = local_data(c)
        record.update(ld)

        # MW data per curve ('gens', 'ngens', 'heights', 'reg',
        # 'torsion_order', 'torsion_structure', 'torsion_gens',
        # 'omega', 'Lvalue', 'sha'):
        record.update(mwdata[i])

        # Galois image data via magma
        galrepdata = get_galrep_data(E, verbose=verbose) # one string with spaces
        record.update(galrepdata)

        data[label] = record

    return data

def make_isogeny_classes(raw_curves, verbose=0, prec=None):
    """raw_curves is a generator yielding short 'raw' curve dicts with fields (as in read_curves_magma()):

    field_label, conductor_norm, conductor_label, conductor_ideal, iso_label, ainvs

    This function computes the isogeny class of each and returns a
    dist whose keys are curve labels, values are curve records.

    """
    data = {}
    for record1 in raw_curves:
        data.update(make_isogeny_class(record1, prec=prec, verbose=verbose))

    return data



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
