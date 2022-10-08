from sys import stdout
import os
from sage.all import Magma, EllipticCurve, NumberField, ZZ, QQ, polygen
from fields import get_field_label, get_IQF_info, get_field_name, ideal_from_IQF_label
from files import read_newform_data, read_missing_levels, BIANCHI_DATA_DIR
from psort import ideal_label, ideal_from_label, primes_iter
from codec import ideal_to_string, old_ideal_label

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
        Elist[i] = EllipticCurve(mag('aInvariants(E);\n').sage()).global_minimal_model(semi_global=True)
    if local_magma:
        mag.quit()
    return Elist

class_number_one_fields = [1, 2, 3, 7, 11, 19, 43, 67, 163]

def magma_search(field, missing_label_file=None, field_info_filename=None, bmf_filename=None, min_norm=None, max_norm=None, outfilename=None, effort=1000, verbose=False):
    r"""
    Uses Magma via EllipticCurveSearch() to search for missing curves (over IQFs given some BMFs).

    INPUT:

    - ``field`` (integer) -- 1, 2, 3, 7, 11;  or 19, 43, 67, 163; or 23, 31, ..., 5, ...

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
    if field in class_number_one_fields:
        if field_info_filename==None:
            field_info_filename = os.path.join(BIANCHI_DATA_DIR, "fieldinfo", "fieldinfo-{}".format(field))
            if verbose:
                print("Using {} for field info".format(field_info_filename))
        K, Plist = get_IQF_info(field_info_filename, 200, verbose)
        field_label = get_field_label(K)
    else:
        x = polygen(QQ)
        if field%4==3:
            K = NumberField(x**2-x+(field+1)//4, 'w')
            field_label = "2.0.{}.1".format(field)
        else:
            K = NumberField(x**2+field, 'w')
            field_label = "2.0.{}.1".format(4*field)
        print("Field {} = {}".format(field_label, K))
        Plist = list(primes_iter(K,maxnorm=ZZ(200)))
    if bmf_filename==None:
        print("Must supply name of a file containing BMFs over {} in {}".format(field, BIANCHI_DATA_DIR))
    else:
        if verbose:
            print("Using {} for newform input".format(bmf_filename))

    if outfilename:
        outfile=open(outfilename, mode="a")
        if verbose:
            print("Using {} for output".format(outfilename))
    output("Field {}\n".format(field_label))
    newforms = read_newform_data(bmf_filename)
    if verbose:
        print("...read newform data finished")
    if missing_label_file==None:
        missing_label_file = bmf_filename
        if verbose:
            print("Using {} for missing labels".format(missing_label_file))

    mag=Magma()
    for level in read_missing_levels(open(missing_label_file)):
        if "." in level:
            N = ideal_from_label(K, level)
        else:
            N = ideal_from_IQF_label(K, level)
        NN = N.norm()
        if min_norm and NN<min_norm:
            continue
        if max_norm and NN>max_norm:
            continue
        goodP = [(i,P) for i,P in enumerate(Plist) if not P.divides(N)]
        level_label = ideal_label(N)
        if verbose:
            print("Missing conductor %s = %s" % (level_label,N))
        nfs = newforms[level]
        for id in nfs.keys():
            nf = nfs[id]
            class_label = "%s-%s" % (level_label,id)
            if verbose:
                print("Working on form %s" % class_label)
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
                print("Found {} curves matching {}: {}".format(len(curves),class_label," ".join([str(E.ainvs()) for E in curves])))
                E = curves[0]
            else:
                print("**********No curve found to match newform {}*************".format(class_label))
                E = None
            # output 3 lines per curve, as expected by the function read_curves_magma:
            output("Conductor {}\n".format(ideal_to_string(N)))
            output("Isogeny_class {}\n".format(class_label))
            if E!=None:
                ainvs = str(list(E.ainvs())).replace("a", "w")
                output("Curve {}\n".format(ainvs))
            else:
                output("No curve found\n")
            if outfilename:
                outfile.flush()

def make_ec_dict(E):
    K = E.base_field()
    N = E.conductor()
    ai = E.ainvs()
    ec = {}
    ec['field_label'] = get_field_label(K)
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
    print("Defining {} prime ideals".format(len(Plist)))
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
    if bmf_filename==None:
        print("Must supply name of a file containing BMFs over {} in {}".format(field, BIANCHI_DATA_DIR))
    else:
        if verbose:
            print("Using {} for newform input".format(bmf_filename))

    if field in class_number_one_fields:
        if field_info_filename==None:
            field_info_filename = os.path.join(BIANCHI_DATA_DIR, "fieldinfo", "fieldinfo-{}".format(field))
        if verbose:
            print("Using {} for field info".format(field_info_filename))
        K, Plist = get_IQF_info(field_info_filename, 200, verbose)
        field_lab = get_field_label(K)
    else:
        field_lab = "2.0.{}.1".format(field)
        x = polygen(QQ)
        if field%4==3:
            K = NumberField(x**2-x+(field+1)//4, 'w')
        else:
            K = NumberField(x**2+field, 'w')
        print("Field {} = {}".format(field_lab, K))
        Plist = list(primes_iter(K,maxnorm=ZZ(200)))
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
        if "." in level:
            N = ideal_from_label(K, level)
        else:
            N = ideal_from_IQF_label(K, level)
        goodP = [(i,P) for i,P in enumerate(Plist) if not P.divides(N)]
        if verbose:
            print("Missing level %s = %s" % (level,N))
            print("{} good primes".format(len(goodP)))
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
            print("aplist has length {}".format(len(aplist)));
            print("nf[ap] has length {}".format(len(nf['ap'])));
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

# functions to create lines of output for files curves.*,
# curve_data.*, isoclass.* from a dictionary containing the relevant
# curve data.

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
        from sage.schemes.elliptic_curves.isogeny_class import IsogenyClass_EC_NumberField
        cl = IsogenyClass_EC_NumberField(Elist[0], reducible_primes=None, algorithm='Billerey', minimal_models=True)
        #cl = Elist[0].isogeny_class()
        perm = dict([(i, cl.index(E)) for i, E in enumerate(Elist)])
        mat = permute_mat(cl.matrix(), perm, True)
        mat = str([list(ri) for ri in mat.rows()]).replace(" ", "")

    output_fields = [ec['field_label'],
                     ec['conductor_label'],
                     ec['iso_label'],
                     str(ec['number']),
                     mat]
    return " ".join(output_fields)

