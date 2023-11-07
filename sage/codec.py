# Functions for coding/decoding data to/from strings

from sage.all import ZZ, QQ, RR, EllipticCurve, prod, KodairaSymbol
from psort import ideal_label
from fields import nf_lookup
from schemas import column_names
from galrep import parse_galrep_data_string

def split_field_label(field_label):
    r"""Return (degree, signature, abs(disc), index) from an LMFDB field label
    """
    d, r, a, i = field_label.split(".")
    d = int(d)
    r = int(r)
    a = ZZ(a)
    i = int(i)
    s = [r, (d-r)//2]
    return (d, s, a, i)

def numerify_iso_label(lab):
    r"""Return the numerical equivalent of an isogeny class letter-code.

    Normally this is the integer>=0 which lab represents in base 26,
    with digits a=0,...,z=25, but if lab starts with 'CM' and the rest
    represents the integer n as above, the numerical version is
    -n-1. (Not -n as then 0 would be ambiguous).  This variant is only
    relevant (currently) over imaginary quadratic fields.  The label
    may be in upper case (currently only over 3.1.23.1), but that
    should not be used over any field for which a CM label is
    possible.
    """
    from sage.databases.cremona import class_to_int
    if 'CM' in lab:
        return -1 - class_to_int(lab[2:])
    else:
        return class_to_int(lab.lower())

# Functions to parse a single line from one of the files curves.*, isoclass.*, mwdata.*, local_data.*, galdata.*
#
# In each case the function returns a full label and a dict whose kets are exactly the relevant table columns
#

# The first 4 columns in curves.*, isoclass.*, mwdata.*, local_data.*
# are the same and define the full label, so we factor out this part.
#


def parse_line_label_cols(L):
    r"""
    Parse the first 4 columns of one line from a curves/isoclass/mwdata/local_data file
    """
    data = L.split()
    record = {}
    record['field_label'] = field_label = data[0]
    degree, signature, abs_disc, _ = split_field_label(field_label)
    record['degree'] = degree
    record['signature'] = signature
    record['abs_disc'] = abs_disc

    conductor_label = data[1]
    iso_label = data[2]
    number = int(data[3])
    short_label = "{}-{}{}".format(conductor_label, iso_label, number)
    label = "{}-{}".format(field_label, short_label)
    short_class_label = "{}-{}".format(conductor_label, iso_label)

    record['conductor_label'] = conductor_label
    record['iso_label'] = iso_label
    record['iso_nlabel'] = numerify_iso_label(iso_label)
    record['number'] = number
    record['short_label'] = short_label
    record['label'] = label
    record['short_class_label'] = short_class_label
    record['class_label'] = "{}-{}".format(field_label, short_class_label)
    return label, record

def parse_curves_line(L):
    r"""Parse one line from a curves file

    14 Columns (see schemas.py):
    'field_label', 'conductor_label', 'iso_label', 'number',
    'conductor_ideal', 'conductor_norm', 'ainvs', 'jinv', 'disc',
    'normdisc', 'equation', 'cm', 'base_change', 'q_curve'
    """
    data = L.split()
    if len(data) != len(column_names['curves']):
        print("curves line {} does not have 12 fields, skipping".format(L))
        return
    label, record = parse_line_label_cols(L)

    record['conductor_ideal'] = data[4]
    record['conductor_norm'] = N = ZZ(data[5])
    record['conductor_norm_factors'] = N.support()

    record['ainvs'] = data[6]
    record['jinv'] = data[7]
    record['disc'] = disc = data[8]
    if "." in disc:
        print("Old disc: {}".format(disc))
        disc = "({})".format(ZZ(RR(disc[1:-1])))
        print("New disc: {}".format(disc))
        record['disc'] = disc
    record['normdisc'] = ZZ(data[9])
    from sage.all import sqrt
    record['root_analytic_conductor'] = sqrt(0.00798504020212804*float(N)**(1.0/float(record['degree']))*float(record['abs_disc']))
    #print('root_analytic_conductor = {}'.format(record['root_analytic_conductor']))

    eqn = data[10]
    # the reason for doing the following is for the unique field
    # 2.2.5.1 where the field generator is not a single character such
    # as 'a' or 'i' but is '\phi', and we don't want to have '\phix'
    # in a latex string (and also do not want any whitespace).
    if "{x}" not in eqn:
        eqn = eqn.replace('x', '{x}').replace('y', '{y}')
    record['equation'] = eqn

    record['cm'] = cm = ZZ(data[11]) if data[11] != '?' else '?'
    # The 'cm_type' column  holds +1 for a curve with rational, -1 for
    # potential, 0 if no CM
    if cm:
        if 'CM' in label:
            record['cm_type'] = +1
        else:
            record['cm_type'] = -1
    else:
        record['cm_type'] = 0
    bc = data[12][1:-1]
    record['base_change'] = [str(lab) for lab in bc.split(",")] if bc else []
    record['q_curve'] = (data[13] == '1')
    return label, record

def parse_isoclass_line(L):
    r"""
    Parse one line from an isoclass file
    """
    data = L.split()
    if len(data) != len(column_names['isoclass']):
        print("isoclass line {} does not have 6 fields, skipping".format(L))
        return
    label, record = parse_line_label_cols(L)

    record['isogeny_matrix'] = mat = [[int(a) for a in r.split(",")]
                                      for r in data[4][2:-2].split("],[")]
    record['class_size'] = len(mat)
    record['class_deg'] = max(max(r) for r in mat)
    record['all_iso_degs'] = dict([[n+1, sorted(list(set(row)))] for n, row in enumerate(mat)]) 
    record['trace_hash'] = ZZ(data[5])

    # NB Every curve in the class has the same 'isogeny_matrix',
    # 'class_size', 'class_deg', and the for the i'th curve in the
    # class (for i=1,2,3,...) its 'isogeny_degrees' column is
    # all_iso_degs[i].

    return label, record

def parse_local_data_line(L):
    r"""
    Parse one line from a local_data file
    """
    data = L.split()
    ncols = len(data)
    if ncols not in [6, 7]:
        print("local_data line {} does not have 6 or 7 fields, skipping".format(L))
        return
    label, record = parse_line_label_cols(L)

    ldstring = "" if (ncols == 6) else data[4]
    ld, ldx = local_data_from_string(ldstring)
    record['local_data'] = ld
    record.update(ldx) # fields 'bad_primes', 'n_bad_primes', 'semistable', 'potential_good_reduction', 'tamagawa_product'

    # The non_min_p column is a list of strings
    # e.g. ['(g)', '(g1,g2)'] while the string in
    # the file will contain [(g),(g1,g2)].

    # Currently the list has 0 or 1 entries but we do not want to rely
    # on this.

    nmp = data[-2]
    #print(nmp)
    record['non_min_p'] = [] if nmp == '[]' else ["("+P+")" for P in nmp[2:-2].split("),(")]
    #print(record['non_min_p'])
    record['minD'] = data[-1]

    return label, record

def parse_mwdata_line(L):
    r"""
    Parse one line from an mwdata file
    """
    data = L.split()
    # if len(data)!=14:
    #     print("mwdata line {} does not have 14 fields, skipping".format(L))
    #     return
    label, record = parse_line_label_cols(L)

    r = data[4]
    record['rank'] = None if r == '?' else int(r)
    r = data[5]
    record['rank_bounds'] = '?' if r == '?' else [int(rb) for rb in r[1:-1].split(",")]
    r = data[6]
    record['analytic_rank'] = None if r == '?' else int(r)
    record['ngens'] = int(data[7])
    gens = data[8]
    record['gens'] = [] if gens == '[]' else gens.replace("[[[", "[[").replace("]]]", "]]").replace("]],[[", "]];[[").split(";")
    record['heights'] = data[9]
    record['reg'] = data[10]
    record['torsion_order'] = nt = int(data[11])
    ts = data[12]
    record['torsion_structure'] = [] if ts == '[]' else [int(t) for t in ts[1:-1].split(",")]
    record['torsion_primes'] = ZZ(nt).prime_divisors()
    record['torsion_gens'] = decode_points_one2many(data[13])

    record['omega'] = None
    record['Lvalue'] = None
    record['sha'] = None

    return label, record

def parse_new_mwdata_line(L):
    r"""
    Parse one line from an mwdata file (with extra columns omega, lvalue, sha)
    """
    #from sage.rings.real_mpfr import RealNumber
    data = L.split()
    if len(data) != len(column_names['mwdata']):
        print("mwdata line {} has only {} fields, skipping".format(L, len(data)))
        return
    label, record = parse_line_label_cols(L)

    def decode_col(col, decoder): # use for columns which may have '?'
        return None if col in ['?', 'None'] else decoder(col)

    record['rank'] = decode_col(data[4], int)
    record['rank_bounds'] = decode_col(data[5], decode_int_list)
    record['analytic_rank'] = decode_col(data[6], int)
    record['ngens'] = int(data[7])
    record['gens'] = decode_points_one2many(data[8])
    record['heights'] = data[9]
    #record['reg'] = decode_col(data[10], RealNumber) if record['ngens'] else 1
    record['reg'] = data[10] if record['ngens'] else 1
    record['torsion_order'] = nt = int(data[11])
    record['torsion_primes'] = ZZ(nt).prime_divisors()
    record['torsion_structure'] = decode_int_list(data[12])
    record['torsion_gens'] = decode_points_one2many(data[13])
    if len(data) == 17:
        #record['omega'] = decode_col(data[14], RealNumber)
        #record['Lvalue'] = decode_col(data[15], RealNumber)
        record['sha'] = decode_col(data[16], int)
        record['omega'] = data[14]
        record['Lvalue'] = data[15]
    else:
        record['omega'] = None
        record['Lvalue'] = None
        record['sha'] = None
    return label, record

def parse_galrep_line(L):
    r"""
    Parse one line from a galrep file
    """
    data = L.split(maxsplit=1)
    record = parse_galrep_data_string("" if len(data) == 1 else data[1])
    record['label'] = label = data[0]
    return label, record

def NFelt(a):
    r""" Returns an NFelt string encoding the element a (in a number field
    K).  This consists of d strings representing the rational
    coefficients of a (with respect to the power basis), separated by
    commas, with no spaces.

    For example the element (3+4*w)/2 in Q(w) gives '3/2,2'.

    If a is already a string, do nothing!
    """
    return a if type(a)==type('') else  ",".join([str(c) for c in list(a)])

def ideal_from_string(K, s, IQF_format=False):
    r"""Returns the ideal of K defined by the string s.  If IQF_format is
    True, this is "[N,c,d]" with N,c,d as in a label, while otherwise
    it is of the form "[N,a,alpha]" where N is the norm, a the least
    positive integer in the ideal and alpha a second generator so that
    the ideal is (a,alpha).  alpha is a polynomial in the variable w
    which represents the generator of K (but may actially be an
    integer).
    """
    #print("s = {}".format(s))
    N, a, alpha = s[1:-1].split(",")
    N = ZZ(N)
    a = ZZ(a)
    if IQF_format:
        d = ZZ(alpha)
        I = K.ideal(N//d, K([a, d]))
    else:
        # 'w' is used for the generator name for all fields for
        # numbers stored in the database
        alpha = alpha.replace('w', str(K.gen()))
        I = K.ideal(a, K(alpha))
    if I.norm() == N:
        return I
    else:
        return "wrong" ## caller must check

def ideal_to_string(I):
    from nfscripts import simplified_gens
    gens = tuple(simplified_gens(I))
    return str(gens).replace(",)", ")").replace(" ", "").replace(str(I.number_field().gen()), 'w')

def ideal_from_IQF_label(K, lab):
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
        lab = lab[1:-1].replace(",", ".")
    a, c, d = [ZZ(x) for x in lab.split(".")]

    a /= d
    P = K.ideal([a, c+d*K.gen()])
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
    return curve_from_strings(K, c['ainvs'])

def parse_NFelt(K, s):
    r"""
    Returns an element of K defined by the string s.
    """
    return K([QQ(c) for c in s.split(",")])

def parse_point(K, s):
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
    substrings separated by semicolons, each substring being a
    comma-separated list of strings representing rational numbers
    representing the NF element with respect to its (power) basis.

    Do nothing if this is already a string.
    """
    return ainvs if type(ainvs)==type('') else ";".join([NFelt(ai) for ai in ainvs])

def ainvs_from_string(K, ainvs):
    r"""Reverse of the previous function: converts a string, representing a
    list of NF elements joined by ";", to a list of actual NF
    elements in K.
    """
    return [parse_NFelt(K, ai) for ai in ainvs.split(";")]

def curve_from_string(K, ainvs):
    r""" Given a number field K and a string, representing a list of 5
    elements, converts these to elements of K and returns the
    elliptic curve with these a-invariants.
    """
    return EllipticCurve(ainvs_from_string(K, ainvs))

# The next two are the old version of the previous two, when we used 5
# separate strings for the a-invariants instead of a single string
# joined by ";" as used in the database itself.

def ainvs_from_strings(K, ainv_string_list):
    r"""
    Reverse of the previous function: converts a list of strings,
    each representing an NF element, to a list of actual NF
    elements in K.
    """
    return [parse_NFelt(K, ai) for ai in ainv_string_list]

def curve_from_strings(K, ainv_string_list):
    r"""
    Given a number field K and a list of 5 strings, each
    representing an NF element, converts these to elements of K
    and returns the elliptic curve with these a-invariants.
    """
    return EllipticCurve(ainvs_from_strings(K, ainv_string_list))

#####################################################################
#
# utility for making a look-up table for converting labels over IQFs
#
#####################################################################

the_labels = {}
field_labels = ['2.0.{}.1'.format(d) for d in [4, 8, 3, 7, 11]]

def convert_ideal_label(K, lab):
    """An ideal label of the form N.c.d is converted to N.i.  Here N.c.d
    defines the ideal I with Z-basis [a, c+d*w] where w is the standard
    generator of K, N=N(I) and a=N/d.  The standard label is N.i where I is the i'th ideal of norm N in the standard ordering.

    NB Only intended for use in coverting IQF labels!  To get the standard label from any ideal I just use ideal_label(I).
    """
    if K in the_labels:
        if lab in the_labels[K]:
            return the_labels[K][lab]
    else:
        the_labels[K] = {}

    comps = lab.split(".")
    # test for labels which do not need any conversion
    if len(comps) == 2:
        return lab
    assert len(comps) == 3
    N, c, d = [int(x) for x in comps]
    a = N//d
    I = K.ideal(a, c+d*K.gen())
    newlab = ideal_label(I)
    #print("Ideal label converted from {} to {} over {}".format(lab, newlab, K))
    the_labels[K][lab] = newlab
    return newlab

def convert_conductor_label(field_label, label):
    """If the field is imaginary quadratic, calls convert_ideal_label, otherwise just return label unchanged.
    """
    if field_label.split(".")[:2] != ['2', '0']:
        return label
    K = nf_lookup(field_label)
    return convert_ideal_label(K, label)

######################################################

def encode_point(P):
    r"""Encodes a point on an elliptic curve over a field of degree d as a
    string representing a 3-list of d-lists of rationals.  Do nothing
    if already a string.
    """
    return P if type(P)==type('') else str([list(c) for c in P]).replace(" ", "")

def encode_points(Plist):
    r"""Converts a list of points into a string encoding a list of 3-lists
    of d-lists of rationals.
    """
    return '[' + ','.join([encode_point(P) for P in Plist]) + ']'

def decode_points_one2many(gens):
    return [] if gens == '[]' else gens.replace("[[[", "[[").replace("]]]", "]]").replace("]],[[", "]];[[").split(";")

def encode_points_many2one(gens):
    return ("["+",".join(gens)+"]").replace(" ", "")

def encode_int_list(L):
    """
    From a list of ints (or reals) return same as a string with no spaces.

    e.g. from [1, 2, 3] return '[1,2,3]'
    """
    return str(L).replace(" ", "")

def decode_int_list(L):
    """
    From a string with no spaces representing a list of ints return the list of ints.

    e.g. from  '[1,2,3]' return [1, 2, 3], from '[]' return []
    """
    return [] if L == '[]' else [int(a) for a in L[1:-1].split(",")]

##########################################################

def local_data_to_string_one_prime(ldp):
    # we do not just join ldp.values() since we want to fix the order
    ldstr = ":".join([str(ldp[k]) for k in ['p', 'normp', 'ord_cond', 'ord_disc', 'ord_den_j', 'red', 'rootno', 'kod', 'cp']])
    ldstr = ldstr.replace(" ", "")
    return ldstr

def local_data_to_string(ld):
    return ";".join([local_data_to_string_one_prime(ldp) for ldp in ld])

#from sage.all import latex
def numerify_kodaira(kod):
    #kod1 = kod
    kod = kod.replace("_", "")
    kod = kod.replace("{", "")
    kod = kod.replace("}", "")
    kod = kod.replace("^", "")
    nkod = KodairaSymbol(kod)._pari_code()
    # kod2 = latex(KodairaSymbol(nkod))
    # if not kod1.replace("{", "").replace("}", "")==kod2.replace("{", "").replace("}", ""):
    #     print("{} --> {} --> {} --> {}".format(kod1, kod, nkod, kod2))
    return nkod

def local_data_from_string_one_prime(s):
    dat = s.split(":")
    return {'p': dat[0], # string
            'normp': int(dat[1]),
            'ord_cond': int(dat[2]),
            'ord_disc': int(dat[3]),
            'ord_den_j': int(dat[4]),
            'red': None if dat[5] == 'None' else int(dat[5]),
            'rootno': '?' if dat[6] == '?' else int(dat[6]),
            'kod': int(dat[7]),
            'cp': int(dat[8])}

def local_data_from_string(s):
    if s:
        ld = [local_data_from_string_one_prime(si) for si in s.split(";")]
    else:
        ld = []
    # ld_extra holds anything else which is per curve not per prime
    ld_extra = {}
    ld_extra['bad_primes'] = bad_primes = [ldp['p'] for ldp in ld if ldp['ord_cond']]
    ld_extra['n_bad_primes'] = len(bad_primes)
    ld_extra['semistable'] = all(ldp['ord_cond'] < 2 for ldp in ld)
    ld_extra['potential_good_reduction'] = all(ldp['ord_den_j'] == 0 for ldp in ld)
    ld_extra['tamagawa_product'] = prod([ldp['cp'] for ldp in ld], 1)
    return ld, ld_extra

##########################################################

num_encoder = str
bool_encoder = lambda b: str(int(b))
string_list_encoder = lambda L: "[" + ",".join(L) + "]"
string_encoder = lambda r: str(r).replace(" ", "") if r != None else '?'
rank_encoder = lambda r: str(r) if r != None else '?'
gal_im_encoder = " ".join

encoders = {'number': num_encoder,
            'conductor_norm': num_encoder,
            'cm': num_encoder,
            'q_curve': bool_encoder,
            'local_data': local_data_to_string,
            'bad_primes': string_list_encoder,
            'non_min_p': string_list_encoder,
            'base_change': string_list_encoder,
            'isogeny_matrix': encode_int_list,
            'trace_hash': num_encoder,
            'rank': rank_encoder,
            'analytic_rank': rank_encoder,
            'sha': rank_encoder,
            'rank_bounds': encode_int_list,
            'gens': encode_points,
            'ngens': num_encoder,
            'torsion_structure': encode_int_list,
            'torsion_gens': encode_points,
            'torsion_order': num_encoder,
            'omega': num_encoder,
            'Lvalue': rank_encoder,
            'reg': rank_encoder,
            'normdisc': num_encoder,
            'ainvs': ainvs_to_string,
            'jinv': NFelt,
            'heights': encode_int_list,
            'galois_images': gal_im_encoder,
            # don't use the next line, it messes up rational coefficients!
            #            'equation': lambda x: x.replace("{","").replace("}","")
}

def get_encoder(col):
    return encoders.get(col, lambda x: x)

def file_line(ftype, c):
    r"""Given a dict containing the data for one curve,
    return the string for one line of a <ftype>.* file.
    """
    if ftype in column_names:
        if ftype == 'mwdata':
            # this fix should not be necessary
            if c['rank'] is None and c['analytic_rank'] is None:
                c['reg'] = None
            # for k in column_names[ftype]:
            #     print("{}: {} --> {} ({})".format(k, c[k], get_encoder(k)(c[k]), type(get_encoder(k)(c[k]))))
        return " ".join([get_encoder(k)(c[k]) for k in column_names[ftype]])
    else:
        raise ValueError("{} is not a valid file type".format(ftype))

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
    torsion_order (torsion order)
    torsion_structure (list of <=2 ints>1, defining the torsion structure)
    torsion_gens (list of points)

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
    hts = str(mwdata['heights']).replace(" ", "")
    reg = str(mwdata['reg'])
    ntors = str(mwdata['torsion_order'])
    torstruct = str(mwdata['torsion_structure']).replace(" ", "")
    tgens = encode_points(mwdata['torsion_gens'])
    output_fields = [c['field_label'], c['conductor_label'], c['iso_label'], str(c['number']),
                     r, rbds, ar, ngens, gens, hts, reg,
                     ntors, torstruct, tgens]
    return " ".join(output_fields)

def make_mwdata_lines(cl):
    r""" return a string for all lines of a mwdata file for one isogeny class.
    """
    return "\n".join([make_mwdata_line(
        {
            'field_label': cl['field_label'],
            'conductor_label': cl['conductor_label'],
            'iso_label': cl['iso_label'],
            'number': i+1,
            'mwdata': mw,
            }
    ) for i, mw in enumerate(cl['mwdata'])])

