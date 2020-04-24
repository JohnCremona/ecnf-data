# coding=utf-8
#
# Functions written to help filling in gaps in the Bianchi (imaginary
# quadratic field) curve tables, using Bianchi newforms data and Magma
# scripts.
#
from sys import stdout
from os import getenv
from sage.all import (polygen, ZZ, QQ, Magma, magma, latex, EllipticCurve, primes,
                      flatten, Primes, legendre_symbol, prod, RR, PowerSeriesRing, O, Integer, srange, sign)
from fields import add_field, field_data, field_label, cm_j_dict, get_IQF_info, get_field_name
from files import read_newform_data, read_missing_levels
from psort import nf_key, primes_of_degree_iter, ideal_from_label
from codec import curve_from_strings, ideal_to_string, old_ideal_label

HOME = getenv("HOME")
bianchi_data_dir = HOME + "/bianchi-data"

# copy of function in lmfdb/scripts/ecnf/hmf_check_find.py
#
def EllipticCurveSearch(K, Plist, N, aplist, effort=1000, mag=None):
    r""" Call Magma's own EllipticCurveSearch() function to find and
    elliptic curve E defined over K with conductor N and ap as in the
    list.

    INPUT:

    - `K` (number field) -- the base field
    - `Plist` (list) -- a list of primes of K
    - `N` (ideal) -- an integral ideal of K
    - `aplist` (list) -- a list of integers a(P) indexed by the P in Plist not dividing N

    OUTPUT:

    A list (possibly empty) of elliptic curves defined over K with
    conductor N and traces of Frobenius given by the a(P).
    """
    # Create a new magma instance for each search:
    local_magma = int(mag is None)
    if local_magma:
        mag = Magma()
    # Define the number field in Magma and the list of primes
    mag.eval("SetGRH();\n")
    mag.eval('SetVerbose("ECSearch",3);\n')
    mag.eval("Qx<x> := PolynomialRing(RationalField());\n")
    name = K.gen()
    pol = K.defining_polynomial()
    mag.eval("Qx<x> := PolynomialRing(RationalField());\n")
    mag.eval("K<%s> := NumberField(%s);\n" % (name, pol))
    mag.eval("OK := Integers(K);\n")
    mag.eval("Plist := [];\n")
    for P in Plist:
        Pgens = P.gens_reduced()
        Pmagma = "(%s)*OK" % Pgens[0]
        if len(Pgens) > 1:
            Pmagma += "+(%s)*OK" % Pgens[1]
        mag.eval("Append(~Plist,%s);\n" % Pmagma)

    mag.eval('SetColumns(0);\n')
    mag.eval('effort := %s;\n' % effort)

    Ngens = N.gens_reduced()
    Nmagma = "(%s)*OK" % Ngens[0]
    if len(Ngens) > 1:
        Nmagma += "+(%s)*OK" % Ngens[1]
    mag.eval('N := %s;' % Nmagma)
    mag.eval('aplist := %s;' % aplist)
    mag.eval('goodP := [P: P in Plist | Valuation(N,P) eq 0];\n')
    mag.eval('goodP := [goodP[i]: i in [1..#(aplist)]];\n')
    try:
        mag.eval('curves := EllipticCurveSearch(N,effort : Primes:=goodP, Traces:=aplist);\n')
        mag.eval('curves := [E: E in curves | &and[TraceOfFrobenius(E,goodP[i]) eq aplist[i] : i in [1..#(aplist)]]];\n')
        mag.eval('ncurves := #curves;')
        ncurves = mag('ncurves;').sage()
    except RuntimeError as arg:
        print("RuntimError in Magma: {}".format(arg))
        if local_magma:
            mag.quit()
        return []
    if ncurves==0:
        if local_magma:
            mag.quit()
        return []
    Elist = [0 for i in range(ncurves)]
    for i in range(ncurves):
        mag.eval('E := curves[%s];\n' % (i+1))
        Elist[i] = EllipticCurve(mag('aInvariants(E);\n').sage())
    if local_magma:
        mag.quit()
    return Elist

def magma_search(field, missing_label_file=None, field_info_filename=None, bmf_filename=None, min_norm=None, max_norm=None, outfilename=None, effort=1000, verbose=False):
    r"""
    Uses Magma via EllipticCurveSearch() to search for missing curves (over IQFs given some BMFs).

    INPUT:

    - ``field`` (integer) -- 1, 2, 3, 7 or 11.

    - ``missing_label_file`` (string) -- filename of file containing
      labels of missing isogeny classes.  If absent, assumes all
      newforms in the newforms file are to be treated.

    - ``field_info_filename`` (string) -- filename of file containing
      field information.  Defaults to
      HOME + "/bianchi-data/fieldinfo/findinfo-%s" % field

    - ``bmf_filename`` (string) -- filename of file containing
      newforms.

    - ``outfilename`` (string, default ``None``) -- name of output file

    - ``effort`` (int, default 1000) -- parameter to Magma's search function

    - ``verbose`` (boolean, default ``False``) -- verbosity flag.

    OUTPUT:

    (To file and/or screen, nothing is returned): Magma commands to
    search for curves given their conductors and Traces of Frobenius,
    the conductors taken from the given file and the traces from the
    associated Bianchi newform file.  The output from the Magma run
    can be parsed by the following function.
    """
    def output(L):
        if outfilename:
            outfile.write(L)
        if verbose:
            stdout.write(L)
    if field_info_filename==None:
        field_info_filename = "%s/fieldinfo/fieldinfo-%s" % (bianchi_data_dir,str(field))
        if verbose:
            print("Using {} for field info".format(field_info_filename))
    if bmf_filename==None:
        print("Must supply name of a file containing BMFs over {} in {}".format(field, bianchi_data_dir))
    else:
        if verbose:
            print("Using {} for newform input".format(bmf_filename))

    K, Plist = get_IQF_info(field_info_filename, 200, verbose)
    field_lab = field_label(K)
    if outfilename:
        outfile=open(outfilename, mode="a")
        if verbose:
            print("Using {} for output".format(outfile))
    newforms = read_newform_data(bmf_filename)
    if verbose:
        print("...read newform data finished")
    if missing_label_file==None:
        missing_label_file = bmf_filename
        if verbose:
            print("Using {} for missing labels".format(missing_label_file))

    bad_labels = []#["16900.0.130-b","16900.0.130-c"]
    mag=Magma()
    for level in read_missing_levels(open(missing_label_file)):
        N = ideal_from_label(K, level)
        NN = N.norm()
        if min_norm and NN<min_norm:
            continue
        if max_norm and NN>max_norm:
            continue
        goodP = [(i,P) for i,P in enumerate(Plist) if not P.divides(N)]
        if verbose:
            print("Missing level %s = %s" % (level,N))
        nfs = newforms[level]
        for id in nfs.keys():
            nf = nfs[id]
            label = "%s-%s" % (level,id)
            if label in bad_labels:
                print("\nIgnoring form %s" % label)
                continue
            else:
                if verbose:
                    print("\nWorking on form %s" % label)
            # Create the array of traces for good primes:
            aplist = [nf['ap'][i] for i,P in goodP if i<len(nf['ap'])]
            # Do the search:
            try:
                curves = EllipticCurveSearch(K, Plist, N, aplist, effort, mag)
            except RuntimeError:
                # Magma throws a run-time error if it finds no curves
                # with the correct traces
                curves = []
            if curves:
                print("Found {} curves matching {}: {}".format(len(curves),label," ".join([str(E.ainvs()) for E in curves])))
                E = curves[0]
            else:
                print("**********No curve found to match newform {}*************".format(label))
                E = None
            if E!=None:
                ec = {}
                ec['field_label'] = field_lab
                ec['conductor_label'] = level
                ec['iso_label'] = id
                ec['number'] = int(1)
                conductor_ideal = E.conductor()
                ec['conductor_ideal'] = ideal_to_string(conductor_ideal,False)
                ec['conductor_norm'] = NN
                ai = E.ainvs()
                ec['ainvs'] = [[str(c) for c in list(a)] for a in ai]
                ec['cm'] = '?'
                ec['base_change'] = []
                output(make_curves_line(ec) + "\n")
                if outfilename:
                    outfile.flush()

def make_ec_dict(E):
    K = E.base_field()
    N = E.conductor()
    ai = E.ainvs()
    ec = {}
    ec['field_label'] = field_label(K)
    ec['conductor_label'] = old_ideal_label(N)
    ec['iso_label'] = 'a'  # placeholder only
    ec['number'] = int(1) # placeholder only
    ec['conductor_ideal'] = ideal_to_string(N,False)
    ec['conductor_norm'] = N.norm()
    ec['ainvs'] = [[str(c) for c in list(a)] for a in ai]
    ec['cm'] = '?'
    ec['base_change'] = []
    return ec


######################################################################
# The functions output_magma_field(), magma_search_script() and
# parse_magma_output() are only useful on a machine which does not
# have Magma installed.  The first two create a magma input file,
# which when run finds curves as requested; the last parses the output
# of such a run.
#
# If Magma is available on the same machine it is easier to use
# EllipticCurveSearch() above.
######################################################################

def output_magma_field(field_label, K, Plist, outfilename=None, verbose=False):
    r"""
    Writes Magma code to a file to define a number field and list of primes.

    INPUT:

    - ``field_label`` (str) -- a number field label

    - ``K`` -- a number field.

    - ``Plist`` -- a list of prime ideals of `K`.

    - ``outfilename`` (string, default ``None``) -- name of file for output.

    - ``verbose`` (boolean, default ``False``) -- verbosity flag.  If
      True, all output written to stdout.

    NOTE:

    Does not assumes the primes are principal.

    OUTPUT:

    (To file and/or screen, nothing is returned): Magma commands to
    define the field `K` and the list `Plist` of primes.
    """
    if outfilename:
        outfile = open(outfilename, mode="w")
    name = K.gen()
    pol = K.defining_polynomial()

    def output(L):
        if outfilename:
            outfile.write(L)
        if verbose:
            stdout.write(L)
    output('print "Field %s";\n' % field_label)
    output("Qx<x> := PolynomialRing(RationalField());\n")
    output("K<%s> := NumberField(%s);\n" % (name, pol))
    output("OK := Integers(K);\n")
    output("Plist := [];\n")
    for P in Plist:
        Pgens = P.gens_reduced()
        Pmagma = "(%s)*OK" % Pgens[0]
        if len(Pgens) > 1:
            Pmagma += "+(%s)*OK" % Pgens[1]
        output("Append(~Plist,%s);\n" % Pmagma)
        # output("Append(~Plist,(%s)*OK);\n" % P.gens_reduced()[0])
    output('effort := 400;\n')
    # output definition of search function:
    output('ECSearch := procedure(class_label, N, aplist);\n')
    output('print "Isogeny class ", class_label;\n')
    output('goodP := [P: P in Plist | Valuation(N,P) eq 0];\n')
    output('goodP := [goodP[i]: i in [1..#(aplist)]];\n')
    output('curves := EllipticCurveSearch(N,effort : Primes:=goodP, Traces:=aplist);\n')
    output('curves := [E: E in curves | &and[TraceOfFrobenius(E,goodP[i]) eq aplist[i] : i in [1..#(aplist)]]];\n')
    output('if #curves eq 0 then print "No curve found"; end if;\n')
    output('for E in curves do;\n ')
    output('a1,a2,a3,a4,a6:=Explode(aInvariants(E));\n ')
    output('printf "Curve [%o,%o,%o,%o,%o]\\n",a1,a2,a3,a4,a6;\n ')
    output('end for;\n')
    output('end procedure;\n')
    output('SetColumns(0);\n')
    if outfilename:
        output("\n")
        outfile.close()

def magma_search_script(field, missing_label_file=None, field_info_filename=None, bmf_filename=None, outfilename=None, verbose=False):
    r"""
    Creates Magma script to search for missing curves.

    INPUT:

    - ``field`` (integer) -- 1, 2, 3, 7 or 11.

    - ``missing_label_file`` (string) -- filename of file containing
      labels of missing isogeny classes.  If absent, assumes all
      newforms in the newforms file are to be treated.

    - ``field_info_filename`` (string) -- filename of file containing
      field information.  Defaults to
      HOME + "/bianchi-data/fieldinfo/findinfo-%s" % field

    - ``bmf_filename`` (string) -- filename of file containing
      newforms.

    - ``outfilename`` (string, default ``None``) -- name of output file

    - ``verbose`` (boolean, default ``False``) -- verbosity flag.

    OUTPUT:

    (To file and/or screen, nothing is returned): Magma commands to
    search for curves given their conductors and Traces of Frobenius,
    the conductors taken from the given file and the traces from the
    associated Bianchi newform file.  The output from the Magma run
    can be parsed by the following function.
    """
    def output(L):
        if outfilename:
            outfile.write(L)
        if verbose:
            stdout.write(L)
    if field_info_filename==None:
        field_info_filename = "%s/fieldinfo/fieldinfo-%s" % (bianchi_data_dir,str(field))
    if bmf_filename==None:
        print("Must supply name of a file containing BMFs over {} in {}".format(field, bianchi_data_dir))
    else:
        if verbose:
            print("Using {} for newform input".format(bmf_filename))

    K, Plist = get_IQF_info(field_info_filename, 200, verbose)
    field_lab = field_label(K)
    if outfilename:
        output_magma_field(field_lab,K,Plist,outfilename)
        if verbose:
            print("...output definition of field and primes finished")
    if outfilename:
        outfile=open(outfilename, mode="a")
    newforms = read_newform_data(bmf_filename)
    if verbose:
        print("...read newform data finished")
    effort = 400;
    output("effort := %s;\n" % effort);
    if missing_label_file==None:
        missing_label_file = bmf_filename

    for level in read_missing_levels(open(missing_label_file)):
        N = ideal_from_label(K, level)
        goodP = [(i,P) for i,P in enumerate(Plist) if not P.divides(N)]
        if verbose:
            print("Missing level %s = %s" % (level,N))
        nfs = newforms[level]
        for id in nfs.keys():
            nf = nfs[id]
            label = "%s.%s" % (level,id)
            if verbose:
                print("\nWorking on form %s" % label)
            # Create the ideal in Magma:
            output('\nprint "Isogeny class %s";\n' % label)
            output("N := (%s)*OK;\n" % N.gens_reduced()[0])
            # Create the array of good primes in Magma:
            output("goodP := [P: P in Plist | Valuation(N,P) eq 0];\n")
            # Create the array of traces for good primes in Magma:
            aplist = [nf['ap'][i] for i,P in goodP if i<len(nf['ap'])]
            output("aplist := %s;\n" % str(aplist))
            output("goodP := [goodP[i]: i in [1..#(aplist)]];\n")
            # Do the search in Magma:
            output("curves := EllipticCurveSearch(N,effort : Primes:=goodP, Traces:=aplist);\n")
            output('if #curves gt 0 then print "Curve found!"; else print "No curve found"; end if;\n')
            output("for E in curves do print aInvariants(E); end for;\n")
    if outfilename:
        outfile.close()

def parse_magma_output(d, infilename, outfilename=None, verbose=False):
    r"""
    Convert Magma search output to a Sage file which can be run to
    provide an iterator through the curves found.

    INPUT:

    - ``d`` (integer) --  1, 2, 3, 7, or 11

    - ``infilename`` (string) -- name of file containing Magma output

    - ``outfilename`` (string, default ``None``) -- name of output file

    - ``verbose`` (boolean, default ``False``) -- verbosity flag.
    """
    disc = [0,-4,-8,-3,0,0,0,-7,0,0,0,-11][d]
    name = get_field_name(disc)
    infile = open(infilename)
    if outfilename:
        outfile=open(outfilename, mode="w")

    def output(L):
        if outfilename:
            outfile.write(L)
        if verbose:
            stdout.write(L)

    if d<3:
        output("K.<%s> = QuadraticField(%s)\n" % (name,-d))
    else:
        output("K.<%s> = NumberField(x^2-x+%s)\n" % (name,(d+1)//4))
    output("extracurves=[\n")
    first=True

    while True:
        try:
            L = infile.next()
        except StopIteration:
            output("]\n\n")
            output("def getcurves():\n")
            output("    for ai in extracurves:\n")
            output("        yield EllipticCurve(K,ai)\n")
            if outfilename:
                outfile.close()
            return

        if L[0]=='[':
            if first:
                first=False
            else:
                output(",\n")
            a1 = infile.next().strip("\n").replace(" ","")
            a2 = infile.next().strip("\n").replace(" ","")
            a3 = infile.next().strip("\n").replace(" ","")
            a4 = infile.next().strip("\n").replace(" ","")
            a6 = infile.next().strip("\n").replace(" ","")
            OL = "[" + "".join([a1,a2,a3,a4,a6]) + "]"
            output(OL)

######################################################################

# 3 functions taken from lmfdb/lmfdb/ecnf/import_ecnf_data.py to
# create lines of output for files curves.*, curve_data.*, isoclass.*
# from a dictionary containing the relevant curve data.

def make_curves_line(ec):
    r""" for ec a curve object from the database, create a line of text to
    match the corresponding raw input line from a curves file.

    ec should be a dict with these fields:

    'field_label', 'conductor_label', 'iso_label', 'number',
    'conductor_ideal', 'conductor_norm', 'ainvs', 'cm', 'base_change'

    Output line fields (13):

    field_label conductor_label iso_label number conductor_ideal conductor_norm a1 a2 a3 a4 a6 cm base_change

    Sample output line:

    2.0.4.1 65.18.1 a 1 [65,18,1] 65 1,1 1,1 0,1 -1,1 -1,0 0 0
    """
    output_fields = [ec['field_label'],
                     ec['conductor_label'],
                     ec['iso_label'],
                     str(ec['number']),
                     ec['conductor_ideal'],
                     str(ec['conductor_norm'])
                     ] + [",".join(t) for t in ec['ainvs']
                          ] + [str(ec['cm']), str(int(len(ec['base_change']) > 0))]
    return " ".join(output_fields)


def make_curve_data_line(ec):
    r""" for ec a curve object from the database, create a line of text to
    match the corresponding raw input line from a curve_data file.

    ec should be a dict with these fields:

    'field_label', 'conductor_label', 'iso_label', 'number'

    and optionally some or all of:

    'rank', 'rank_bounds', 'analytic_rank', 'gens', 'sha_an'

    Output line fields (9+n where n is the 8th); all but the first 4
    are optional and if not known should contain"?" except that the 8th
    should contain 0.

    field_label conductor_label iso_label number rank rank_bounds analytic_rank ngens gen_1 ... gen_n sha_an

    Sample output line:

    2.0.4.1 65.18.1 a 1 0 ? 0 0 ?
    """
    rk = '?'
    if 'rank' in ec:
        rk = str(ec['rank'])
    rk_bds = '?'
    if 'rank_bounds' in ec:
        rk_bds = str(ec['rank_bounds']).replace(" ", "")
    an_rk = '?'
    if 'analytic_rank' in ec:
        an_rk = str(ec['analytic_rank'])
    ngens = '0'
    gens_str = []
    if 'gens' in ec:
        gens_str = ["[" + ":".join([c for c in P]) + "]" for P in ec['gens']]
        ngens = str(len(gens_str))
    sha = '?'
    if 'sha_an' in ec:
        sha = str(int(ec['sha_an']))

    output_fields = [ec['field_label'],
                     ec['conductor_label'],
                     ec['iso_label'],
                     str(ec['number']),
                     rk, rk_bds, an_rk,
                     ngens] + gens_str + [sha]
    return " ".join(output_fields)


nfcurves = None # keep pyflakes happy: the function below might need
                # it to be set to the nfcurves collection.
def make_isoclass_line(ec):
    r""" for ec a curve object from the database, create a line of text to
    match the corresponding raw input line from an isoclass file.

    ec should be a dict with these fields:

    'field_label', 'conductor_label', 'iso_label', 'number'

    and optionally:

    'isogeny_matrix'

    Output line fields (15):

    field_label conductor_label iso_label number isogeny_matrix

    Sample output line:

    2.0.4.1 65.18.1 a 1 [[1,6,3,18,9,2],[6,1,2,3,6,3],[3,2,1,6,3,6],[18,3,6,1,2,9],[9,6,3,2,1,18],[2,3,6,9,18,1]]
    """
    mat = ''
    if 'isogeny_matrix' in ec:
        mat = str(ec['isogeny_matrix']).replace(' ', '')
    else:
        print("Making isogeny matrix for class %s" % ec['label'])
        from lmfdb.ecnf.isog_class import permute_mat
        from lmfdb.ecnf.WebEllipticCurve import FIELD
        K = FIELD(ec['field_label'])
        curves = nfcurves.find({'field_label': ec['field_label'],
                                'conductor_label': ec['conductor_label'],
                                'iso_label': ec['iso_label']}).sort('number')
        Elist = [EllipticCurve([K.parse_NFelt(x) for x in c['ainvs']]) for c in curves]
        cl = Elist[0].isogeny_class()
        perm = dict([(i, cl.index(E)) for i, E in enumerate(Elist)])
        mat = permute_mat(cl.matrix(), perm, True)
        mat = str([list(ri) for ri in mat.rows()]).replace(" ", "")

    output_fields = [ec['field_label'],
                     ec['conductor_label'],
                     ec['iso_label'],
                     str(ec['number']),
                     mat]
    return " ".join(output_fields)

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

    add_field(K, field_label=field_label)
    Kdata = field_data[K]

    # CM curves are Q-curves:
    if cm_j_dict[jE]:
        if verbose:
            print("Yes: E is CM")
        return True

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
        aP0 = ap(E,Plist[0])
        for P in Plist[1:]:
            aP = ap(E,P)
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
    if Kdata['is_galois']:
        autos = Kdata['autos']
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
                print("Yes: class contains all conjugate j-invariants")
            else:
                print("No: class does not contain all conjugate j-invariants")
            return t

    return '?'


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
    if any([ld.bad_reduction_type()==0 for ld in E.local_data()]):
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
                   'kod': str(latex(ld.kodaira_symbol())),
                   'cp': int(ld.tamagawa_number())}
                  for ld in Eld]
    return E_local_data, E.non_minimal_primes(), E.minimal_discriminant_ideal()
