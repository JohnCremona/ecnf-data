# -*- coding: utf-8 -*-
r"""Functions to find ranks and generators of elliptic curves by
 interfacing with Magma

This version reads curves from a curves.* file and outputs a curvedata.* file
"""
from __future__ import print_function

from sage.all import EllipticCurve
from fields import nf_lookup
#from curves import curve_from_strings

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
    K = E.base_field()

    def convert_point(P):
        return E([K(c.sage()) for c in P.Eltseq()])
    if verbose:
        print("calling magma...")
    from sage.all import magma
    if HeightBound is None:
        MWSI = magma(E).MordellWeilShaInformation(nvals=3)
    else:
        MWSI = magma(E).MordellWeilShaInformation(nvals=3, HeightBound=HeightBound)
    if verbose:
        print("...done.")
    rank_bounds = MWSI[0].sage()
    gens = [convert_point(P) for P in MWSI[1]]
    sha_bounds = dict(MWSI[2].sage())
    if gens and test_saturation:
        if verbose:
            print("testing that Magma's generators are saturated...")
        newgens, index, newreg = E.saturation(gens, verbose)
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
    for nstep in range(ncurves):  # upper bound for number if iterations needed
        for i in range(ncurves):
            for j in range(ncurves):
                if Qlists[i] != [] and (maps[i][j] != 0) and Qlists[j] == []:
                    # print("Mapping from %s to %s at step %s" % (i,j,nstep))
                    phi = maps[i][j]
                    p = phi.degree()  # a prime
                    Qlists[j] = [maps[i][j](P) for P in Qlists[i]]
                    # now do p-saturation (if possible)
                    try:
                        E = Qlists[j][0].curve()
                        pts, index, reg = E.saturation(Qlists[j], one_prime=p)
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

    A list of pairs [rank_bounds, gens], one for each curve in the class.
    """
    # source = find_source(Cl.isogenies())
    # adiscs = [E.discriminant().norm().abs() for E in Cl.curves]
    # print("Abs disc list: %s" % adiscs)
    ss = [len(str(E.ainvs())) for E in Cl.curves]
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
    from sage.all import union
    n = len(maps)
    mat = [[int((maps[i][j] != 0) or (i == j)) for j in range(n)] for i in range(n)]
    # print("mat = %s" % mat)
    children = [[j for j in range(n) if mat[i][j]] for i in range(n)]
    for nstep in range(n):  # upper bound for number if iterations needed
        for i in range(n):
            for j in children[i]:
                children[i] = union(children[i], children[j])
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

    A list of pairs [rank_bounds, gens], one for each class.
    """
    Cl = curves[0].isogeny_class()
    MWI = MWInfo_class(Cl, HeightBound=HeightBound, test_saturation=test_saturation, verbose=verbose)
    # Now we must map the points to the correct curves!

    n = len(Cl.curves)
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

def encode_point(P):
    r"""
    Converts a point into a string encoding a 3-list of d-lists of rationals
    """
    return str([list(c) for c in list(P)]).replace(" ","")

def encode_points(Plist):
    r"""
    Converts a list of points into a list of lists of 3 lists of d lists of strings
    """
    return [encode_point(P) for P in Plist]

def get_generators(iso_class, test_saturation=False, verbose=False):
    r"""Given an isogeny class, finds the rank (or bounds) and
    generators.

    iso_class is a dict as yielded by the read_classes() function, with keys:
    field_label, N_label, N_def, iso_label, curves
    of which the first 4 are strings and the last is a list of EllipticCurves.

    Returns iso_class with an extra key 'mwdata' whose value is a list
    of data, for each curve in the class, the data being a dict with
    keys 'rank_bounds', 'rank', 'gens', 'heights', 'reg'.  If the rank
    bounds are not equal then key 'rank' is missing.

    """
    class_label = "-".join([iso_class['field_label'], iso_class['N_label'], iso_class['iso_label']])
    Es = iso_class['curves']
    if verbose:
        print("Curves in class %s: %s" % (class_label, [E.ainvs() for E in Es]))
    mwi = MWInfo_curves(Es, HeightBound=2, test_saturation=test_saturation, verbose=verbose)
    if verbose:
        print("MW data: %s" % mwi)

    iso_class['mwdata'] = []
    for E, mw in zip(Es, mwi):
        data = {}
        data['rank_bounds'] = [int(r) for r in mw[0]]
        if mw[0][0] == mw[0][1]:
            data['rank'] = int(mw[0][0])
        gens = mw[1]
        data['gens'] = encode_points(gens)
        data['heights'] = [str(P.height()) for P in gens]
        data['reg'] = E.regulator_of_points(gens)
        if verbose:
            print("MW data for E={}:\n{}".format(E.ainvs(),data))
        iso_class['mwdata'].append(data)

    return iso_class

def make_curve_data_line(c):
    r""" return a string for one line of  a curve_data file.

    c is a dict with keys:
    field_label, conductor_label, iso_label, number, mwdata
    where mwdata is a dict with keys
    rank (int or '?')
    rank_bounds (list of 2 ints)
    analytic_rank (int or '?')
    ngens (int: 0 means we have no gens, whatever the rank)
    gens (list of points)
    sha_an (real)

    Output line fields (9+n where n is the 8th); all but the first 4
    are optional and if not known should contain"?" except that the 8th
    should contain 0.

    field_label conductor_label iso_label number rank rank_bounds analytic_rank ngens gen_1 ... gen_n sha_an

    Sample output line:

    2.0.4.1 2053.1809.1 a 1 2 [2,2] ? 2 [[0,0],[-1,0],[1,0]] [[2,0],[2,0],[1,0]] ?
    """
    mwdata = c['mwdata']
    gens = mwdata.get('gens',[])
    ngens = str(len(gens))
    r = str(mwdata['rank']) if 'rank' in mwdata else '?'
    rbds = str(mwdata['rank_bounds']).replace(" ", "")
    ar = str(mwdata['analytic_rank']) if 'analytic_rank' in mwdata else '?'
    sha = str(int(mwdata['sha_an'])) if 'sha_an' in mwdata else '?'

    cond_lab = c['N_label']
    output_fields = [c['field_label'], c['N_label'], c['iso_label'], str(c['number']),
                     r, rbds, ar, ngens] + gens + [sha]
    return " ".join(output_fields)

def make_curve_data_lines(cl):
    r""" return a string for all lines of  a curve_data file for one isogeny class.
    """
    return "\n".join([make_curve_data_line(
        {
            'field_label': cl['field_label'],
            'N_label': cl['N_label'],
            'iso_label': cl['iso_label'],
            'number': i+1,
            'mwdata': mw,
            }
    ) for i,mw in enumerate(cl['mwdata'])])

def make_curve_data(curves_filename, curve_data_filename, min_cond_norm=None, max_cond_norm=None,
                    test_saturation=False, verbose=False):
    r""" Retrieves curves from a curves file
    with conductor norm between given bounds (optional), finds their
    ranks (or bounds) and generators, and outputs a curve_data file.
    """
    cd_out = open(curve_data_filename, 'w')
    for cl in read_classes(curves_filename):
        NN = cl['N_norm']
        if min_cond_norm and NN<min_cond_norm:
            if verbose:
                print("Skipping class as conductor norm < {}".format(min_cond_norm))
            continue
        if max_cond_norm and NN>max_cond_norm:
            if verbose:
                print("Skipping rest of file as conductor norms >= {} > {}".format(NN,max_cond_norm))
            #continue
            break

        if True:#verbose:
            class_label = "-".join([cl['field_label'],cl['N_label'],cl['iso_label']])
            print("Processing class {}".format(class_label))
        cl = get_generators(cl, test_saturation=test_saturation, verbose=verbose)
        cdlines = make_curve_data_lines(cl)
        if verbose:
            print(cdlines)
        cd_out.write(cdlines+"\n")
    cd_out.close()
