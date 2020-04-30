# -*- coding: utf-8 -*-
r"""Functions to find ranks and generators of elliptic curves by
 interfacing with Magma

This version reads curves from a curves.* file and outputs a mwdata.* file

Format of mwdata file supercedes that of curve_data file, including
also heights and regulator and torsion data.

1. field_label
2. conductor_label
3. iso_label
4. number
5: rank                   int or ?
6: rank_bounds      [int,int] or ?
7: analytic_rank          int or ?
8: ngens (number of generators of infinite order) int
9: gens  (list of ngens lists of 3 lists of d rationals)
10: heights         list of ngens reals
11: regulator       real
12: ntorsion         int
13: torstruct       [int,int]
14: torgens (list of 0,1 or 2 lists of 3 lists of d rationals)

"""
from __future__ import print_function
from sage.all import magma, union
from files import read_classes, read_classes_new, parse_mwdata_line
from codec import encode_points

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
        return E([K([ci.sage() for ci in c.Eltseq()]) for c in P.Eltseq()])
    if verbose:
        print("calling magma on E={} over {}...".format(E.ainvs(),K))
    mE = magma(E)
    if HeightBound is None:
        MWSI = mE.MordellWeilShaInformation(nvals=3)
    else:
        MWSI = mE.MordellWeilShaInformation(nvals=3, HeightBound=HeightBound)
    ar = int(mE.AnalyticRank())
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
            'analytic_rank': ar,
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
                        try:
                            pts, index, reg = E.saturation(Qlists[j], one_prime=p, debug=True)
                        except RuntimeError:
                            print("RuntimeError!")
                            print("Base field K = {}".format(E.base_field()))
                            print("E = {}".format(E.ainvs()))
                            print("Points = {}".format(Qlists[j]))
                            print("one_prime = {}".format(p))
                            index=1
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
    return [[MWI['rank_bounds'], pts, MWI['analytic_rank']] for pts in map_points(Cl.isogenies(), source, MWI['gens'], verbose)]


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

    A list of triples [rank_bounds, gens, analytic_rank], one for each class.
    """
    Cl = curves[0].isogeny_class()
    MWI = MWInfo_class(Cl, HeightBound=HeightBound, test_saturation=test_saturation, verbose=verbose)
    # Now we must map the points to the correct curves!

    n = len(Cl.curves)
    if not n==len(curves):
        print("MWInfo_curves list of curves had length {} but isogeny class from 1st curve has length {}".format(len(curves),n))
    fixed_MWI = [None for _ in range(n)]  # just to set length
    for i in range(n):
        E = curves[i]
        j = Cl.index(E)  # checks for isomorphism, not just equality
        iso = Cl.curves[j].isomorphism_to(E)
        # print("(i,j)=(%s,%s)" % (i,j))
        fixed_MWI[i] = [MWI[j][0], [iso(P) for P in MWI[j][1]], MWI[j][2]]

    # Check we have it right:
    assert all([all([P in curves[i] for P in fixed_MWI[i][1]]) for i in range(n)])
    return fixed_MWI

def get_generators(iso_class, test_saturation=False, verbose=False):
    r"""Given an isogeny class, finds the rank (or bounds) and
    generators and torsion data for each curve in the class.

    iso_class is a dict as yielded by the read_classes() function, with keys:
    field_label, N_label, N_def, iso_label, curves
    of which the first 4 are strings and the last is a list of EllipticCurves.

    Returns iso_class with an extra key 'mwdata' whose value is a list
    of data, for each curve in the class, the data being a dict with
    keys

    'rank_bounds', 'rank', 'gens', 'hts', 'reg', 'ntors', 'torstruct', 'tgens'.

    If the rank bounds are not equal then key 'rank' is missing.

    """
    class_label = "-".join([iso_class['field_label'], iso_class['N_label'], iso_class['iso_label']])
    Es = iso_class['curves']
    if verbose:
        print("Curves in class %s: %s" % (class_label, [E.ainvs() for E in Es]))
    try:
        mwi = MWInfo_curves(Es, HeightBound=2, test_saturation=test_saturation, verbose=verbose)
        if verbose:
            print("MW data: %s" % mwi)
    except RuntimeError as e:
        # We can still try for the anayltic rank:
        ar = int(magma(Es[0]).AnalyticRank())
        print(e)
        print("Unable to compute rank bounds for {}, but analytic rank = {}".format(class_label, ar))
        mwi = [[None, [], ar] for E in Es]

    iso_class['mwdata'] = []
    for E, mw in zip(Es, mwi):
        data = {}
        if mw[0]: # else something went wrong and we have no rank bounds (or any gens)
            data['rank_bounds'] = [int(r) for r in mw[0]]
            if mw[0][0] == mw[0][1]:
                data['rank'] = int(mw[0][0])
        data['gens'] = gens = mw[1]
        data['hts'] = [P.height() for P in gens]
        data['reg'] = E.regulator_of_points(gens) if gens else 1
        data['analytic_rank'] = mw[2]

        # get torsion order, structure and generators:
        torgroup = E.torsion_subgroup()
        data['ntors'] = torgroup.order()
        data['torstruct'] = list(torgroup.invariants())
        data['tgens'] = [P.element() for P in torgroup.gens()]

        if verbose:
            print("MW data for E={}:\n{}".format(E.ainvs(),data))
        iso_class['mwdata'].append(data)

    return iso_class

def make_mwdata_line(c):
    r"""return a string for one line of a mwdata file.

    c is a dict with keys:
    field_label, conductor_label, iso_label, number, mwdata
    where mwdata is a dict with keys
    rank (int or '?')
    rank_bounds (list of 2 ints)
    analytic_rank (int or '?')
    ngens (int: 0 means we have no gens, whatever the rank)
    gens (list of points)
    ntors (torsion order)
    torstruct (list of <=2 ints>1, defining the torsion structure)
    tgens (list of points)

    Output line fields (14); all but the first 4 are optional. If
    rank, rank_bounds or analytic_rank are not known they must be "?".

    field_label conductor_label iso_label number rank rank_bounds
    analytic_rank ngens gens heights regulator ntors torstruct torgens

    Sample output line:

    2.0.4.1 2053.1809.1 a 1 2 [2,2] ? 2 [[[0,0],[-1,0],[1,0]];[[2,0],[2,0],[1,0]]] 1 [] []

    """
    mwdata = c['mwdata']
    r = str(mwdata['rank']) if 'rank' in mwdata else '?'
    rbds = str(mwdata['rank_bounds']).replace(" ", "") if 'rank_bounds' in mwdata else '?'
    ar = str(mwdata['analytic_rank']) if 'analytic_rank' in mwdata else '?'
    ngens = str(len(mwdata['gens']))
    gens = encode_points(mwdata['gens'])
    hts = str(mwdata['hts']).replace(" ", "")
    reg = str(mwdata['reg'])
    ntors = str(mwdata['ntors'])
    torstruct = str(mwdata['torstruct']).replace(" ", "")
    tgens = encode_points(mwdata['tgens'])
    output_fields = [c['field_label'], c['N_label'], c['iso_label'], str(c['number']),
                     r, rbds, ar, ngens, gens, hts, reg,
                     ntors, torstruct, tgens]
    return " ".join(output_fields)

def make_mwdata_lines(cl):
    r""" return a string for all lines of a mwdata file for one isogeny class.
    """
    return "\n".join([make_mwdata_line(
        {
            'field_label': cl['field_label'],
            'N_label': cl['N_label'],
            'iso_label': cl['iso_label'],
            'number': i+1,
            'mwdata': mw,
            }
    ) for i,mw in enumerate(cl['mwdata'])])

def make_mwdata(curves_filename, mwdata_filename, label=None,
                min_cond_norm=None, max_cond_norm=None,
                test_saturation=False, verbose=False):
    r"""Retrieves curves from a curves file
    with conductor norm between given bounds (optional), finds their
    ranks (or bounds) and generators, and outputs an mwdata file.

    If label is given, it should be a short isogeny class label, and
    then only that class will be run.  Otherwise, the minimum and
    maximum conductor may optionally be given.  Otherwise all the
    curves (isogeny classes) in the in put file are processed.

    """
    with open(mwdata_filename, 'w', 1) as mw_out:
        for cl in read_classes_new(curves_filename):
            short_class_label = "-".join([cl['N_label'],cl['iso_label']])
            class_label = "-".join([cl['field_label'],cl['N_label'],cl['iso_label']])
            if label:
                if label!=short_class_label:
                    if verbose:
                        print("Skipping {}".format(short_class_label))
                    continue
            NN = cl['N_norm']
            if min_cond_norm and NN<min_cond_norm:
                # if verbose:
                #     print("Skipping class as conductor norm < {}".format(min_cond_norm))
                continue
            if max_cond_norm and NN>max_cond_norm:
                if verbose:
                    print("Skipping rest of file as conductor norms >= {} > {}".format(NN,max_cond_norm))
                break

            print("Processing class {}".format(class_label))
            try:
                cl = get_generators(cl, test_saturation=test_saturation, verbose=verbose)
                mwlines = make_mwdata_lines(cl)
                if verbose:
                    print(mwlines)
                mw_out.write(mwlines+"\n")
            except RuntimeError as e:
                print("caught RuntimeError: {}".format(e))

def add_analytic_ranks(curves_filename, mwdata_filename, suffix='x', verbose=False):
    r"""Retrieves curves from a curves file and mwdata from the mwdata
     file.  Computes analytic ranks and rewrites the mwdata file
     adding the suffix to its filename.

    This is a one-off since the orginal mwdata file code forgot to
    compute and output analytic ranks.
    """
    ar_table = {}
    n = 0
    for cl in read_classes_new(curves_filename):
        short_class_label = "-".join([cl['N_label'],cl['iso_label']])
        class_label = "-".join([cl['field_label'],cl['N_label'],cl['iso_label']])
        ar = int(magma(cl['curves'][0]).AnalyticRank())
        ar_table[class_label] = ar
        if verbose:
            print("Processing class {}: analytic rank = {}".format(class_label, ar))
        n += 1
    print("Finished computing analytic ranks for {} classes".format(n))
    #print(ar_table)
    with open(mwdata_filename) as mw_in, open(mwdata_filename+suffix, 'w', 1) as mw_out:
        for L in mw_in.readlines():
            label, record = parse_mwdata_line(L)
            if record['analytic_rank'] is None:
                class_label = record['class_label']
                ar = ar_table[class_label]
                if verbose:
                    print("Updating analytic rank of {} to {}".format(class_label,ar))
                data = L.split()
                data[6] = str(ar)
                L = " ".join(data)
                if verbose:
                    print("New mwdata line: {}".format(L))
                mw_out.write(L + "\n")
            else:
                mw_out.write(L)
