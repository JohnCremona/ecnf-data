# coding=utf-8
#
# Functions written to help filling in gaps in the Bianchi (imaginary
# quadratic field) curve tables, using Bianchi newforms data and Magma
# scripts.
#
import sys
from sage.all import polygen, ZZ, QQ, NumberField, PolynomialRing, Magma, EllipticCurve

field_names=dict([(-4,'i'),(-8,'t'),(-3,'w')])
def get_field_name(disc):
    r"""
    Returns the name of the generator for the field of given discriminant.

    INPUT:

    - ``disc`` (integer)-- a field discriminant

    OUTPUT:

    'w', 'i', 't' for discriminants -3.-4.-8, else 'a'.
    """
    return field_names.get(disc,'a')

def field_data(s):
    r"""
    Returns full field data from field label.

    INPUT:

    - ``s`` (string)-- an LMFDB field label

    OUTPUT:

    List containing the input string, degree, signature (list), absolute discriminant.
    """
    deg, r1, abs_disc, n = [int(c) for c in s.split(".")]
    sig = [r1, (deg-r1)//2]
    return [s, deg, sig, abs_disc]

def field_from_label(lab):
    r"""
    Returns a number field from its LMFDB label.

    INPUT:

    - ``s`` (string)-- an LMFDB field label

    OUTPUT:

    A number field.
    """
    dummy, deg, sig, abs_disc = field_data(lab)
    x = polygen(QQ)
    if deg==2:
        d = ZZ(abs_disc)
        if sig[0]==0: d=-d
        t = d%4
        assert t in [0,1]
        pol = x**2 - t*x + (t-d)/4
    elif lab=='3.1.23.1':
        pol = x**3 - x**2 +1
    else:
        raise NotImplementedError("cannot yet handle field %s" % lab)
    K = NumberField(pol, 'a')
    print "Created field from label %s: %s" % (lab,K)
    return K

def field_label(K):
    r"""
    Returns the LMFDB label of a number field.

    *** Only works when the label's last component is 1 ***

    INPUT:

    - ``K`` -- a number field

    OUTPUT:

    (string) the LMFDB label of K.
    """
    d = K.degree()
    r = K.signature()[0] # number of real embeddings
    D = K.discriminant().abs()
    return "{}.{}.{}.1".format(d,r,D)
    return "%s.%s.%s.1" % (d,r,D)

def ideal_from_label(K,lab):
    r"""
    Returns an ideal in quadratic field K from its label.

    INPUT:

    - ``K`` -- a quadratic number field

    - ``lab`` (string) -- label of an ideal in K

    OUTPUT:

    The ideal defined by the label.  Labels have the form '[N,c,d]'
    where `N` is the norm of the ideal and the HNF of the ideal is
    `\left<a,c+d\alpha\right>` with `a=N/d` and `\alpha` the standard
    integral generator of `K`.
    """
    if '[' in lab:
        lab = lab[1:-1].replace(",",".")
    a,c,d = [ZZ(x) for x in lab.split(".")]

    a /= d
    P = K.ideal([a,c+d*K.gen()])
    return P


# HNF of an ideal I in a quadratic field

def ideal_HNF(I):
    r"""
    Returns an HNF triple defining the ideal I in a quadratic field
    with integral basis [1,w].

    This is a list [a,b,d] such that [a,c+d*w] is a Z-basis of I, with
    a,d>0; c>=0; N = a*d = Norm(I); d|a and d|c; 0 <=c < a.
    """
    N = I.norm()
    (a, c), (b, d) = [[ZZ(x) for x in row] for row in I.pari_hnf().python()]
    assert a > 0 and d > 0 and N == a * d and d.divides(a) and d.divides(b) and 0 <= c < a
    return [a, c, d]

# Label of an ideal I in a quadratic field: string formed from the
# Norm and HNF of the ideal

def old_ideal_label(I):
    r"""
    Returns the HNF-based label of an ideal I in a quadratic field
    with integral basis [1,w].  This is the string 'N.c.d' where
    [a,c,d] is the HNF form of I and N=a*d=Norm(I).
    """
    a, c, d = ideal_HNF(I)
    return "%s.%s.%s" % (a * d, c, d)

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

def get_field_info(field_info_filename, maxpnorm=200, verbose=False):
    r"""
    Returns a number field and ordered list of primes.

    INPUT:

    - ``field_info_filename`` (string) -- name of data file.

    - ``maxpnorm`` (integer) -- bound on norm of primes to return.

    - ``verbose`` (boolean, default False) -- verbosity flag.

    OUTPUT:

    Tuple of a number field and a list of prime ideals, ordered as in the data file.
    """
    field_info_file = file(field_info_filename)
    Plist=[]
    for L in field_info_file.readlines():
        if "Q" in L: # first line
            Zx = PolynomialRing(ZZ,'x')
            poly = Zx(L.split()[-1][:-1])
            name = get_field_name(poly.discriminant())
            K = NumberField(poly,name)
            if verbose:
                print("Field is %s" % K)
        else:
            nm,lab,gen,p,e,deg = L.split()
            nm = ZZ(nm)
            if nm>maxpnorm:
                break
            p = ZZ(p)
            e = ZZ(e)
            deg = ZZ(deg)
            P = ideal_from_label(K,lab)
            assert P.norm()==nm
            #print("Prime %s with char %s, degree %s, label %s, norm %s" % (P,p,deg,lab,nm))
            Plist += [P]
    field_info_file.close()
    return K, Plist

def nf_filename_from_D(absD):
    d=absD if absD%4 else absD//4
    return "/home/jec/bianchi-data/nflist/nflist.%s.1-20000" % d

def read_newform_data(nf_filename, verbose=False):
    r"""
    Returns a dict containing Bianchi newforms read from a file.

    INPUT:

    - ``nf_filename`` (string) -- name of file containing Bianchi
      newform data.  Either starts with "nflist" or "newforms"; the
      latter has additional fields which we ignore here.

    - ``verbose`` (boolean, default False) -- verbosity flag.

    OUTPUT:

    dict with keys level labels (strings), values dicts with keys
    newforms labels (strings 'a', 'b', etc), values dicts with keys
    'gen_str' (string representing generator of the level as a
    principal ideal), 'sign' (sign of functional equation of the
    L-function), 'aq' (list of Atkin-Lehner eigenvalues), 'ap' (list
    of L-function coefficients of prime ideals in standard order).

    """
    nf_file = file(nf_filename)
    old_fmt =  "nflist" in nf_filename
    print("file has {} format".format('old' if old_fmt else 'new'))
    newforms = {}
    for L in nf_file.readlines():
        #if verbose:
        #    print("raw input: %s" % L)
        if old_fmt:
            label, gen, sfe, loverp, ALs, aplist = L.split()
            level, letter = label.split("-")
        else:
            field_label, level, letter, gen, wt, bc, cm, sfe, loverp, ALs, poly, aplist = L.split()

        if not level in newforms:
            newforms[level] = {}
        nfs = newforms[level]
        nfs[letter] = nf = {}
        nf['gen_str'] = gen[1:-1]
        nf['sign'] = int(sfe)
        nf['aq'] = [int(e) for e in ALs[1:-1].split(",")]
        nf['ap'] = [int(e) for e in aplist[1:-1].split(",")]
        if verbose:
            print("newform data: %s" % L)

    nf_file.close()
    return newforms

def read_missing_levels(infile):
    r"""
    Yields level labels from a file of newform/isogeny class labels, without repeats.
    """
    levels = []
    for L in infile.readlines():
        if "nflist" in infile:
            level = L.split("-")[0]
        else:
            level = L.split()[1]
        if not level in levels:
            levels += [level]
            yield level

bianchi_data_dir = "/home/jec/bianchi-data"

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
    except RuntimeError, arg:
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

def magma_search(field, missing_label_file=None, field_info_filename=None, nf_filename=None, min_norm=None, max_norm=None, outfilename=None, effort=1000, verbose=False):
    r"""
    Uses Magma via EllipticCurveSearch() to search for missing curves.

    INPUT:

    - ``field`` (integer) -- 1, 2, 3, 7 or 11.

    - ``missing_label_file`` (string) -- filename of file containing
      labels of missing isogeny classes.  If absent, assumes all
      newforms in the newforms file are to be treated.

    - ``field_info_filename`` (string) -- filename of file containing
      field information.  Defaults to
      "/home/jec/bianchi-data/fieldinfo/findinfo-%s" % field

    - ``nf_filename`` (string) -- filename of file containing
      newforms.  Defaults to
      "/home/jec/bianchi-data/nflist.%s.1-20000" % field

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
            sys.stdout.write(L)
    if field_info_filename==None:
        field_info_filename = "%s/fieldinfo/fieldinfo-%s" % (bianchi_data_dir,str(field))
        if verbose:
            print("Using {} for field info".format(field_info_filename))
    if nf_filename==None:
        nf_filename = "%s/nflist/nflist.%s.1-20000" % (bianchi_data_dir,str(field))
        if verbose:
            print("Using {} for newform input".format(nf_filename))

    K, Plist = get_field_info(field_info_filename, 200, verbose)
    field_lab = field_label(K)
    if outfilename:
        outfile=file(outfilename, mode="a")
        if verbose:
            print("Using {} for output".format(outfile))
    newforms = read_newform_data(nf_filename)
    if verbose:
        print("...read newform data finished")
    if missing_label_file==None:
        missing_label_file = nf_filename
        if verbose:
            print("Using {} for missing labels".format(missing_label_file))

    bad_labels = []#["16900.0.130-b","16900.0.130-c"]
    mag=Magma()
    for level in read_missing_levels(file(missing_label_file)):
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
        outfile = file(outfilename, mode="w")
    name = K.gen()
    pol = K.defining_polynomial()

    def output(L):
        if outfilename:
            outfile.write(L)
        if verbose:
            sys.stdout.write(L)
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

def magma_search_script(field, missing_label_file=None, field_info_filename=None, nf_filename=None, outfilename=None, verbose=False):
    r"""
    Creates Magma script to search for missing curves.

    INPUT:

    - ``field`` (integer) -- 1, 2, 3, 7 or 11.

    - ``missing_label_file`` (string) -- filename of file containing
      labels of missing isogeny classes.  If absent, assumes all
      newforms in the newforms file are to be treated.

    - ``field_info_filename`` (string) -- filename of file containing
      field information.  Defaults to
      "/home/jec/bianchi-data/fieldinfo/findinfo-%s" % field

    - ``nf_filename`` (string) -- filename of file containing
      newforms.  Defaults to
      "/home/jec/bianchi-data/nflist.%s.1-10000" % field

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
            sys.stdout.write(L)
    if field_info_filename==None:
        field_info_filename = "%s/fieldinfo/fieldinfo-%s" % (bianchi_data_dir,str(field))
    if nf_filename==None:
        nf_filename = "%s/nflist.%s.1-10000" % (bianchi_data_dir,str(field))

    K, Plist = get_field_info(field_info_filename, 200, verbose)
    if outfilename:
        output_magma_field(K,Plist,outfilename)
        if verbose:
            print("...output definition of field and primes finished")
    if outfilename:
        outfile=file(outfilename, mode="a")
    newforms = read_newform_data(nf_filename)
    if verbose:
        print("...read newform data finished")
    effort = 400;
    output("effort := %s;\n" % effort);
    if missing_label_file==None:
        missing_label_file = nf_filename

    for level in read_missing_levels(file(missing_label_file)):
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
    infile = file(infilename)
    if outfilename:
        outfile=file(outfilename, mode="w")

    def output(L):
        if outfilename:
            outfile.write(L)
        if verbose:
            sys.stdout.write(L)

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

#####################################################################
#
# utility for making a look-up table for converting labels over IQFs
#
#####################################################################

from psort import ideal_label

the_labels = {}
field_labels = ['2.0.{}.1'.format(d) for d in [4,8,3,7,11]]
the_fields = dict([(lab,field_from_label(lab)) for lab in field_labels])
print(the_fields)

def convert_ideal_label(K, lab):
    """An ideal label of the form N.c.d is converted to N.i.  Here N.c.d
    defines the ideal I with Z-basis [a, c+d*w] where w is the standard
    generator of K, N=N(I) and a=N/d.  The standard label is N.i where I is the i'th ideal of norm N in the standard ordering.

    NB Only intended for use in coverting IQF labels!  To get the standard label from any ideal I just use ideal_label(I).
    """
    global the_labels
    if K in the_labels:
        if lab in the_labels[K]:
            return the_labels[K][lab]
    else:
        the_labels[K] = {}

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
    the_labels[K][lab] = newlab
    return newlab

def label_conversion_table(infile, outfile):
    out = file(outfile, mode='w')
    for L in file(bianchi_data_dir + "/ideals/" + infile).readlines():
        field, ideal = L.split()
        label = convert_ideal_label(the_fields[field],ideal)
        out.write(' '.join([field, ideal, label])+'\n')
    out.close()

