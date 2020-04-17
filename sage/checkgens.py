# coding=utf-8

# Given files curves.x and curve_data.x with the same curves, check
# that the generators in the latter lie on the correct curve,
# i.e. that the labels are consistent.
#
# Usage:
#        sage: %runfile checkgens.py
#        sage: check_gens(x)
#
# Will first report how many curves were read from curves.x.  Then,
# for each line in curve_data.x, will report if the label does not
# match any label in the curves file, or if the label matches but the
# points specified in the curve_data file do not lie on the curve
# specified in the curves file with the same label.  Otherwise, no
# news is good news.  No warning is given for curves in the curves
# file which do not appear in the curve_data file (e.g. incomplete
# data in the latter).

from sage.all import QQ, ZZ, polygen, NumberField, EllipticCurve
from fields import nf_lookup

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

def parse_point(K,s):
    r"""
    Returns an list of 3 elements of K defined by the string s.

    Example: K=Q(a) quadratic
             s = '[[-302/9,-16/3],[3098/27,-685/27],[1,0]]'

    returns [-302/9-(16/3)*a, 3098/27-(685/27)*a, 1]
    """
    return [K([QQ(c) for c in coord.split(",")]) for coord in s[2:-2].split('],[')]

fields = {}
# This function only works for quadratic fields, where the label
# defined the field uniquely and easily without having to look up in
# the LMFDB!
def field_from_label(lab):
        if lab in fields:
                return fields[lab]
        dummy, deg, sig, abs_disc = field_data(lab)
        d = ZZ(abs_disc)
        if sig[0]==0: d=-d
        x = polygen(QQ)
        t = d%4
        assert t in [0,1]
        pol = x**2 - t*x + (t-d)/4
        K = NumberField(pol, 'a')
        fields[lab] = K
        print("Created field from label {}: {}".format(lab,K))
        return K

def parse_curves_line(L):
        data = L.split()
        if len(data)!=13:
            print("line {} does not have 13 fields, skipping".format(L))
            return ('', None)
        K = nf_lookup(data[0])
        ainvs = [parse_NFelt(K,ai) for ai in data[6:11]]
        E = EllipticCurve(ainvs)

        field_label = data[0]       # string
        conductor_label = data[1]   # string
        iso_label = data[2]         # string
        number = int(data[3])       # int
        short_label = "%s-%s%s" % (conductor_label, iso_label, str(number))
        label = "%s-%s" % (field_label, short_label)
        return (label, E)

def parse_curve_data_line(L):
        data = L.split()
        ngens = int(data[7])
        if len(data)!=9+ngens:
            print("line {} does not have 9 fields (excluding gens), skipping".format(L))
        field_label = data[0]       # string
        conductor_label = data[1]   # string
        iso_label = data[2]         # string
        number = int(data[3])       # int
        short_label = "%s-%s%s" % (conductor_label, iso_label, str(number))
        label = "%s-%s" % (field_label, short_label)
        return (label, data[8:8+ngens])

def read_curves(infile):
    curves = {}
    for L in open(infile).readlines():
        label, E = parse_curves_line(L)
        if label:
            yield (label, E)
        else:
            print("line {} does not have 13 fields, skipping".format(L))
            continue

def check_gens(suffix, verbose=False):
    curves_file = "curves.%s" % suffix
    curvedata_file = "curve_data.%s" % suffix
    cfile = open(curves_file)
    cdfile = open(curvedata_file)
    all_good = True
    bad_curves = []
    n = 0
    m = 0
    while True:
        L1 = cfile.readline()
        if not L1:
            break
        L2 = cdfile.readline()
        if not L2:
            break
        n += 1
        if verbose and n%100==0:
            print("{} curves read".format(n))
        lab2, pts = parse_curve_data_line(L2)
        if not pts:
            continue
        m += 1
        label, E = parse_curves_line(L1)
        if label!=lab2:
            print("label mismatch! {} from {} but {} from {}".format(label,curves_file, lab2, curvedata_file))
            cfile.close()
            cdfile.close()
            return

        K = E.base_field()
        for pt in pts:
            try:
                assert E(parse_point(K,pt))
                if verbose:
                    print("{} OK on {}".format(pt,label))
            except:
                print("Bad point {} for curve {}".format(pt,label))
                if not label in bad_curves:
                    bad_curves.append(label)
                all_good = False
    print("Processed {} curves of which {} had any points".format(n,m))
    if all_good:
        print("All generators check OK")
    else:
        print("Bad generators for {}".format(bad_curves))
    cfile.close()
    cdfile.close()
