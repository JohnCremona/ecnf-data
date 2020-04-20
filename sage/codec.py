# Functions for coding/decoding data to/from strings

from sage.all import ZZ, QQ, EllipticCurve
from fields import nf_lookup

# The next 2 functions are copied from lmfdb/ecnf/WebEllipticCurve.py

def ideal_from_string(K,s, IQF_format=False):
    r"""Returns the ideal of K defined by the string s.  If IQF_format is
    True, this is "[N,c,d]" with N,c,d as in a label, while otherwise
    it is of the form "[N,a,alpha]" where N is the norm, a the least
    positive integer in the ideal and alpha a second generator so that
    the ideal is (a,alpha).  alpha is a polynomial in the variable w
    which represents the generator of K (but may actially be an
    integer).  """
    #print("ideal_from_string({}) over {}".format(s,K))
    N, a, alpha = s.split(".")
    N = ZZ(N)
    a = ZZ(a)
    if IQF_format:
        d = ZZ(alpha)
        I = K.ideal(N//d, K([a, d]))
    else:
        # 'w' is used for the generator name for all fields for
        # numbers stored in the database
        alpha = alpha.encode().replace('w',str(K.gen()))
        I = K.ideal(a,K(alpha.encode()))
    if I.norm()==N:
        return I
    else:
        return "wrong" ## caller must check

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

def ideal_from_IQF_label(K,lab):
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

def curve_from_data(c):
    K = nf_lookup(c['field_label'])
    return curve_from_strings(K,c['ainvs'])

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

def ainvs_to_string(ainvs):
        r"""
        Convert a list of n NF elements to a string consisting of n
        substrings separated by spaces, each substring being a
        comma-separated list of strings representing rational numbers
        representing the NF element with respect to its (power) basis.
        """
        return " ".join([",".join([str(c) for c in list(ai)]) for ai in ainvs])

def ainvs_from_strings(K, ainv_string_list):
        r"""
        Reverse of the previous function: converts a list of strings,
        each representing an NF element, to a list of actual NF
        elements in K.
        """
        return [parse_NFelt(K,ai) for ai in ainv_string_list]

def curve_from_strings(K, ainv_string_list):
        r"""
        Given a number field K and a list of 5 strings, each
        representing an NF element, converts these to elements of K
        and returns the elliptic curve with these a-invariants.
        """
        return EllipticCurve(ainvs_from_strings(K,ainv_string_list))

#####################################################################
#
# utility for making a look-up table for converting labels over IQFs
#
#####################################################################

from psort import ideal_label

the_labels = {}
field_labels = ['2.0.{}.1'.format(d) for d in [4,8,3,7,11]]

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

def convert_conductor_label(field_label, label):
    """If the field is imaginary quadratic, calls convert_ideal_label, otherwise just return label unchanged.
    """
    if field_label.split(".")[:2] != ['2','0']:
        return label
    K = nf_lookup(field_label)
    return convert_ideal_label(K,label)

######################################################

def encode_point(P):
    r"""
    Encodes a point on an elliptic curve over a field of degree d as a string representing a 3-list of d-lists of rationals
    """
    return str([list(c) for c in P]).replace(" ","")

def encode_points(Plist):
    r"""
    Converts a list of points into a string encoding a list of 3-lists of d-lists of rationals
    """
    return '[' + ','.join([encode_point(P) for P in Plist]) + ']'

