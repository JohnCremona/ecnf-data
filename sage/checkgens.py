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
        print "Created field from label %s: %s" % (lab,K)
        return K

def read_curves(infile):
    curves = {}
    for L in file(infile).readlines():
        data = L.split()
        if len(data)!=13:
            print "line %s does not have 13 fields, skipping" % L
            continue
        K = field_from_label(data[0])
        ainvs = [parse_NFelt(K,ai) for ai in data[6:11]]
        E = EllipticCurve(ainvs)

        field_label = data[0]       # string
        conductor_label = data[1]   # string
        iso_label = data[2]         # string
        number = int(data[3])       # int
        short_label = "%s-%s%s" % (conductor_label, iso_label, str(number))
        label = "%s-%s" % (field_label, short_label)
        curves[label] = E

    return curves

def check_gens(suffix):
    curves_file = "curves.%s" % suffix
    curvedata_file = "curve_data.%s" % suffix
    curves = read_curves(curves_file)
    print "Read %s curves" % len(curves)
    for L in file(curvedata_file).readlines():
        data = L.split()
        ngens = int(data[7])
        if not ngens:
            continue
        if len(data)!=9+ngens:
            print "line %s does not have 9 fields (excluding gens), skipping" % line
        field_label = data[0]       # string
        conductor_label = data[1]   # string
        iso_label = data[2]         # string
        number = int(data[3])       # int
        short_label = "%s-%s%s" % (conductor_label, iso_label, str(number))
        label = "%s-%s" % (field_label, short_label)
        try:
            E = curves[label]
        except KeyError:
            print "%s is not the label of a curve read from the curves file" % label
            continue
        K = E.base_field()
        genstringlists = [g[1:-1].split(":") for g in data[8:8+ngens]]
        try:
            gens = [E([parse_NFelt(K,c) for c in cc]) for cc in genstringlists]
        except:
            print "Bad generators %s for curve %s" % (genstringlists,label)
        #print "Curve %s: all %s gens OK" % (label,ngens)
