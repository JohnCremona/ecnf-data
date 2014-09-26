# coding=utf-8
#
# Functions written to help filling in gaps in the Bianchi (imaginary
# quadratic field) curve tables, using Bianchi newforms data and Magma
# scripts.
#
import sys

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
        pol = x^2 - t*x + (t-d)/4
    elif lab=='3.1.23.1':
        pol = x**3 - x**2 +1
    else:
        raise NotImplementedError("cannot yet handle field %s" % lab)
    K = NumberField(pol, 'a')
    print "Created field from label %s: %s" % (lab,K)
    return K

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
    a,c,d = [ZZ(x) for x in lab[1:-1].split(",")]
    a /= d
    P = K.ideal([a,c+d*K.gen()])
    return P


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

def output_magma_field(K,Plist,outfilename=None, verbose=False):
    r"""
    Writes Magma code to a file to define a number field and list of primes.

    INPUT:

    - ``K`` -- a number field.

    - ``Plist`` -- a list of prime ideals of `K`.

    - ``outfilename`` (string, default ``None``) -- name of file for output.

    - ``verbose`` (boolean, default ``False``) -- verbosity flag.  If
      True, all output written to stdout.

    NOTE:

    Assumes the primes are principal: only the first generator is used
    in the Magma ideal construction.

    OUTPUT:

    (To file and/or screen, nothing is returned): Magma commands to
    define the field `K` and the list `Plist` of primes.
    """
    if outfilename:
        outfile=file(outfilename, mode="w")
    disc = K.discriminant()
    name = get_field_name(disc)
    pol = K.defining_polynomial()
    def output(L):
        if outfilename:
            outfile.write(L)
        if verbose:
            sys.stdout.write(L)
    output("Qx<x> := PolynomialRing(RationalField());\n")
    output("K<%s> := NumberField(%s);\n" % (name,pol))
    output("OK := Integers(K);\n")
    output("Plist := [];\n")
    for P in Plist:
        output("Append(~Plist,(%s)*OK);\n" % P.gens_reduced()[0])
    if outfilename:
        outfile.close()

def nf_filename(absD):
    d=absD if absD%4 else absD//4
    return "/home/jec/bianchi-data/nflist.%s.1-10000" % d

def read_newform_data(nf_filename, verbose=False):
    r"""
    Returns a dict containing Bianchi newforms read from a file.

    INPUT:

    - ``nf_filename`` (string) -- name of file containing Bianchi newform data.

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
    newforms = {}
    for L in nf_file.readlines():
        if verbose:
            print("raw input: %s" % L)
        label, gen, sfe, loverp, ALs, aplist = L.split()
        level, letter = label.split(".")
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
    return newforms

def read_missing_levels(infile):
    r"""
    Yields level labels from a file of newform/isogeny class labels, without repeats.
    """
    levels = []
    for L in infile.readlines():
        level = L.split(".")[0]
        if not level in levels:
            levels += [level]
            yield level

bianchi_data_dir = "/home/jec/bianchi-data"

def magma_search_script(field, missing_label_file, outfilename=None, verbose=False):
    r"""
    Creates Magma script to search for missing curves.

    INPUT:

    - ``field`` (integer) -- 1, 2, 3, 7 or 11.

    - ``missing_label_file`` (string) -- filename of file containing
      labels of missing isogeny classes.

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
    field_info_filename = "%s/fieldinfo-%s" % (bianchi_data_dir,str(field))
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
