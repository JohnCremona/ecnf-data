# -*- coding: utf-8 -*-
r"""Functions to find ranks and generators of elliptic curves by
 interfacing with Magma

"""
from __future__ import print_function
from sage.all import RealField, ZZ
from nfscripts import torsion_data, global_period, analytic_rank_and_lvalue

AR_DEGREE_BOUND = 6 # do not compute analytic ranks over field of degree larger than this

# cache values of  RR(K.discriminant().abs()).sqrt()

# NB previous versons divided this factor by 2**(K.signature()[1]) =
# 2**r2 where r2 = #complex places, but instead we adjust the global
# period by doubling the period at complex places.

field_factors = {}
def get_field_factor(K, RR):
    global field_factors
    if K not in field_factors:
        field_factors[K] = RR(K.discriminant().abs()).sqrt()
    return field_factors[K]

def MWShaInfo(E, HeightBound=None, test_saturation=False, verbose=False):
    r"""
    Interface to Magma's MordellWeilShaInformation function

    INPUT:

    - E : an elliptic curve defined over a number field (including Q)

    OUTPUT:

    a dict with keys rank_bounds, gens, sha_bounds where

    - rank_bounds is a list of 2 integers [r1,r2] such r1 <= rank(E) <= r2

    - gens is a list of r1 independent points on E

    - sha_bounds is a dict with keys positive integers n, values rank bounds for Sha[n].

    EXAMPLE::

        sage: E = EllipticCurve('5077a1')
        sage: MWShaInfo(E)
        {'rank_bounds': [3, 3],
        'gens': [(3 : -4 : 1), (2 : -1 : 1), (-1 : 3 : 1)],
        'sha_bounds': {2: [0, 0]}}

    """
    from magma import get_magma
    K = E.base_field()

    def convert_point(P):
        return E([K([ci.sage() for ci in c.Eltseq()]) for c in P.Eltseq()])
    if verbose:
        print("calling magma on E={} over {}...".format(E.ainvs(), K))
    magma = get_magma()
    mE = magma(E)
    if HeightBound is None:
        MWSI = mE.MordellWeilShaInformation(nvals=3)
    else:
        MWSI = mE.MordellWeilShaInformation(nvals=3, HeightBound=HeightBound)
    if verbose:
        print("...done.")
    rank_bounds = MWSI[0].sage()
    gens = [convert_point(P) for P in MWSI[1]]
    sha_bounds = dict(MWSI[2].sage())
    if gens:
        maxp = 0 if test_saturation else 100
        if verbose:
            print("testing that Magma's generators are saturated", end="")
            if maxp:
                print(f" (at primes up to {maxp} only)...")
            else:
                print("...")
        newgens, index, _ = E.saturation(gens, verbose=0, max_prime=maxp)
        if index > 1:
            # Must print this even if not verbose!
            print("Magma's generators for curve %s were not saturated!  index = %s" % (E.ainvs(), index))
            gens = newgens
        else:
            if verbose:
                print("... and they are!")
    return {'rank_bounds': rank_bounds,
            'gens': gens,
            'sha_bounds': sha_bounds,
           }


def map_points(maps, source, Plist, verbose=False):
    r"""Given a matrix of isogenies and a list of points on one curve (with
    index 'source'), returns a list of their images on each other curve.
    Since the isogenies only exist for some i,j pairs, we need to know
    an index 'source' such that following the maps from that curve will
    cover the class, and the initial points must be on that curve.

    We assume that the original points are saturated; after mapping
    under a p-isogeny the images may not be p-saturated so additional
    p-saturation is done.
    """
    ncurves = len(maps)
    if ncurves == 1:
        return [Plist]
    Qlists = [[]] * ncurves
    Qlists[source] = Plist
    if not Plist:
        return Qlists
    nfill = 1
    # print("Qlists = %s" % Qlists)
    # while True: // OK if input satisfies the conditions, but otherwise would loop for ever
    for _ in range(ncurves):  # upper bound for number of iterations needed
        for i in range(ncurves):
            for j in range(ncurves):
                if Qlists[i] != [] and (maps[i][j] != 0) and Qlists[j] == []:
                    # print("Mapping from %s to %s at step %s" % (i,j,nstep))
                    phi = maps[i][j]
                    p = phi.degree()  # a prime
                    # Fix Sage bug
                    for P in Qlists[i]:
                        if hasattr(P, '_order'):
                            del P.__dict__['_order']
                    Qlists[j] = [maps[i][j](P) for P in Qlists[i]]
                    # now do p-saturation (if possible)
                    try:
                        E = Qlists[j][0].curve()
                        try:
                            pts, index, _ = E.saturation(Qlists[j], one_prime=p, debug=True)
                        except RuntimeError:
                            print("RuntimeError!")
                            print("Base field K = {}".format(E.base_field()))
                            print("E = {}".format(E.ainvs()))
                            print("Points = {}".format(Qlists[j]))
                            print("one_prime = {}".format(p))
                            index = 1
                        if index > 1:
                            Qlists[j] = E.lll_reduce(pts)[0]
                            if verbose:
                                print("%s-saturation needed on curve %s, gaining index %s" % (p, list(E.ainvs()), index))
                        else:
                            if verbose:
                                print("image points on curve %s already %s-saturated" % (j, p))
                    except AttributeError:
                        print("Unable to %s-saturate, use a newer Sage version!" % p)
                    # print("...now Qlists = %s" % Qlists)
                    nfill += 1
                    if nfill == ncurves:
                        return Qlists
    if nfill < ncurves:
        raise RuntimeError("In map_points, failed to cover the class!")


def MWInfo_class(Cl, HeightBound=None, test_saturation=False, verbose=False):
    r"""
    Get MW info for all curves in the class

    INPUT:

    - Cl: an isogeny class

    OUTPUT:

    A list of triples [rank_bounds, gens, analytic_rank], one for each curve in the class.
    """
    # source = find_source(Cl.isogenies())
    # adiscs = [E.discriminant().norm().abs() for E in Cl.curves]
    # print("Abs disc list: %s" % adiscs)
    ss = [len(str(E.ainvs())) for E in Cl.curves]
    # n = len(Cl.curves)
    # source = 1 if n==2 else ss.index(min(ss))
    source = ss.index(min(ss))
    if verbose:
        print("Using curve %s to find points" % list(Cl.curves[source].ainvs()))
    MWI = MWShaInfo(Cl.curves[source], HeightBound=HeightBound, test_saturation=test_saturation, verbose=verbose)
    return [[MWI['rank_bounds'], pts] for pts in map_points(Cl.isogenies(), source, MWI['gens'], verbose)]


def find_source(maps):
    r""" maps is an nxn array representing a directed graph on n vertices,
    with some entries 0 and some non-zero.  There will be at least one
    index i such that starting from vertex i you can reach every
    vertex, and we return such an i.  For example if
    maps=[[0,0],[1,0]] then vertex 1 is good but vertex 0 is not!
    """
    n = len(maps)
    mat = [[int((maps[i][j] != 0) or (i == j)) for j in range(n)] for i in range(n)]
    # print("mat = %s" % mat)
    children = [[j for j in range(n) if mat[i][j]] for i in range(n)]
    for _ in range(n):  # upper bound for number if iterations needed
        for i in range(n):
            for j in children[i]:
                children[i] = list(Set(children[i]).union(Set(children[j])))
                if len(children[i]) == n:
                    # print("source = %s" % i)
                    return i
    raise ValueError("find_source problem with mat=%s: at end, children = %s but no source found!" % (mat, children))


def MWInfo_curves(curves, HeightBound=None, test_saturation=False, verbose=False):
    r""" Get MW info for all curves in the list; this is a list of all
    curves in a class in some order, where we do not have the maps
    between them, so we recompute the class with maps and take into
    account the order of the curves.

    INPUT:

    - Cl: an isogeny class

    OUTPUT:

    A list of triples [rank_bounds, gens, analytic_rank], one for each class.
    """
    from sage.schemes.elliptic_curves.isogeny_class import IsogenyClass_EC_NumberField
    Cl = IsogenyClass_EC_NumberField(curves[0], reducible_primes=None, algorithm='Billerey', minimal_models=True)
    #Cl = curves[0].isogeny_class()
    MWI = MWInfo_class(Cl, HeightBound=HeightBound, test_saturation=test_saturation, verbose=verbose)
    # Now we must map the points to the correct curves!

    n = len(Cl.curves)
    if n != len(curves):
        print("MWInfo_curves list of curves had length {} but isogeny class from 1st curve has length {}".format(len(curves), n))
    fixed_MWI = [None for _ in range(n)]  # just to set length
    for i in range(n):
        E = curves[i]
        j = Cl.index(E)  # checks for isomorphism, not just equality
        iso = Cl.curves[j].isomorphism_to(E)
        # print("(i,j)=(%s,%s)" % (i,j))
        fixed_MWI[i] = [MWI[j][0], [iso(P) for P in MWI[j][1]]]

    # Check we have it right:
    assert all([all([P in curves[i] for P in fixed_MWI[i][1]]) for i in range(n)])
    return fixed_MWI

def compute_mwdata(iso_class, test_saturation=False, backend='Magma', verbose=False, prec=None):
    r"""Given an isogeny class, finds the rank (or bounds) and
    generators and torsion data for each curve in the class.

    backend is 'Magma' or 'pari' and controls who computes analytic rank and L-value.

    iso_class is a dict as yielded by the read_classes() function, with keys:
    field_label, conductor_label, conductor_ideal, iso_label, curves
    of which the first 4 are strings and the last is a list of EllipticCurves.

    Returns a list of dicts, one for each curve in the class, each with keys

    'analytic_rank', 'rank_bounds', 'rank', 'gens', 'heights', 'reg', 'torsion_order', 'torsion_structure', 'torsion_gens'.

    - If the rank bounds are not equal then key 'rank' is missing.
    - The analytic rank and L-value are only computed if the base field has degree at most AR_DEGREE_BOUND.

    """
    if prec is None:
        RR = RealField()
        prec = RR.precision()
    else:
        RR = RealField(prec)
    if verbose:
        print(f"In compute_mwdata() with {prec=}, {backend=}, {RR=}")

    class_label = "-".join([iso_class['field_label'], iso_class['conductor_label'], iso_class['iso_label']])
    Es = iso_class['curves']
    E = Es[0]
    K = E.base_field()
    if verbose:
        print("Curves in class %s: %s" % (class_label, [E.ainvs() for E in Es]))
    mwi = MWInfo_curves(Es, HeightBound=2, test_saturation=test_saturation, verbose=verbose)
    if verbose:
        print("MW data: %s" % mwi)

    # analytic rank (for small degree only):
    if K.degree() <= AR_DEGREE_BOUND:
        ar, lval = analytic_rank_and_lvalue(E, prec=prec, backend=backend, verbose=verbose)
    else:
        ar = lval = None

    mwdata = []
    for E, mw in zip(Es, mwi):
        data = {}
        data['analytic_rank'] = ar # may be None
        data['Lvalue'] = lval      # may be None
        data['rank'] = None
        if mw[0]: # else something went wrong and we have no rank bounds (or any gens)
            data['rank_bounds'] = [int(r) for r in mw[0]]
            if mw[0][0] == mw[0][1]:
                data['rank'] = int(mw[0][0])
        data['gens'] = gens = mw[1]
        data['ngens'] = ngens = len(gens)
        data['heights'] = [P.height(prec) for P in gens]
        # compute regulator if either rank lower and upper bounds
        # agree, or lower rank bound = analytic rank:
        if data['rank'] is not None or (ar is not None and mw[0][0]==ar):
            reg = E.regulator_of_points(gens, prec) if gens else 1
        else:
            reg = None
        data['reg'] = reg

        # get torsion order, structure and generators:
        data.update(torsion_data(E))

        # allow for the scaling in the Neron-Tate height
        # pairing: for BSD we need non-normalised heights and
        # normalization divides every height by K.degree(), so
        # the regulator we need has to be multiplied by
        # K.degree()**rank.
        if reg:
            NTreg = reg * K.absolute_degree()**ngens
            if verbose:
                print("Neron-Tate regulator = {}".format(NTreg))
        else:
            NTreg = None

        # compute omega

        # find scaling factor in case we don't have a global minimal model
        minDnorm = E.minimal_discriminant_ideal().norm()
        modelDnorm = E.discriminant().norm().abs()
        fac = (modelDnorm/minDnorm).nth_root(12) # will be exact
        if fac != 1 and verbose:
            print("Not a global minimal model")
            print("Scaling factor = {}".format(fac))
        data['omega'] = omega = global_period(E, fac, prec=prec)
        if verbose:
            print("omega = {}".format(omega))

        nt = data['torsion_order']
        if verbose:
            print("Torsion order = {}".format(nt))

        tamagawa_product = E.tamagawa_product()
        if verbose:
            print("Tamagawa product = {}".format(tamagawa_product))

        # Compute analytic sha, only if (1) we compted the L-value and
        # analytic rank, and (2) we computed the regulator:
        if lval and NTreg:
            Rsha = lval * nt**2  / (NTreg * tamagawa_product * omega)
            Kfactor = get_field_factor(K, RR)
            if verbose:
                print("Field factor = {}".format(Kfactor))
            Rsha *= Kfactor
            data['sha'] = sha = Rsha.round()
            if verbose:
                print("Approximate analytic Sha = {}, rounds to {}".format(Rsha, sha))
            if sha == 0 or (sha-Rsha).abs() > 0.0001 or not ZZ(sha).is_square():
                if not verbose:
                    print("Approximate analytic Sha = {}, rounds to {}".format(Rsha, sha))
                print("****************************Not good! 0 or non-square or not close to a positive integer!")

        else:
            sha = None
            if verbose:
                print("Unable to compute analytic Sha, since we have not computed the rank or special L-value")
        data['sha'] = sha

        if verbose:
            print("MW data for E={}:\n{}".format(E.ainvs(), data))
        mwdata.append(data)

    return mwdata

def get_generators(iso_class, test_saturation=False, backend='Magma', verbose=False, prec=None):
    r"""Given an isogeny class, finds the rank (or bounds) and
    generators and torsion data for each curve in the class.

    backend is 'Magma' or 'pari' and controls who computes analytic rank and L-value.

    iso_class is a dict as yielded by the read_classes() function, with keys:
    field_label, conductor_label, conductor_ideal, iso_label, curves
    of which the first 4 are strings and the last is a list of EllipticCurves.

    Returns iso_class with an extra key 'mwdata' whose value is a list
    of data, for each curve in the class, the data being a dict with
    keys

    'rank_bounds', 'rank', 'gens', 'heights', 'reg', 'torsion_order', 'torsion_structure', 'torsion_gens'.

    If the rank bounds are not equal then key 'rank' is missing.

    """
    iso_class['mwdata'] = compute_mwdata(iso_class, test_saturation, verbose, prec)
    return iso_class

