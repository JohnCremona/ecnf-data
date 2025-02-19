# We have to look up number fields in the database from their labels,
# but only want to do this once for each label, so we will maintain a
# dict of label:field pairs:

from sage.all import QQ, ZZ, PolynomialRing, NumberField, cm_j_invariants_and_orders
from psort import primes_iter

# nf_table is a dict with a field label as key and value the
# associated number field.  It is created by the function
# read_all_fields() which reads from a text file, default
# 'ecnf_fields' which has one field per line, the label then a space
# then the list of ciefficients of the standard defining polynomial.

# NB to process curves over a new field it is essential that the field
# is assed to the ecnf_fields file.

nf_table = {}

# subfield_table is a dict with a field label as key and value a list
# of (proper, nontrivial) subfield labels, possibly empty. It can be
# created using make_subfields() or read from a file (default
# 'ecnf_subfields').

subfield_table = {}

special_names = {'2.0.4.1': 'i',
                 '2.2.5.1': 'phi',
                 '4.0.125.1': 'zeta5',
                 }

field_names=dict([(-4,'i'),(-8,'t'),(-3,'w')])
def get_field_name(disc):
    r"""
    Returns the name of the generator for the field of given discriminant.

    INPUT:

    - ``disc`` (integer)-- a field discriminant

    OUTPUT:

    'w', 'i', 't' for discriminants -3,-4,-8, else 'a'.
    """
    return field_names.get(disc,'a')

def read_all_fields(ffilename='ecnf_fields'):
    ffile = open(ffilename)
    Qx = PolynomialRing(QQ, 'x')
    for line in ffile.readlines():
        label, coeffs = line.split()
        poly = Qx([ZZ(c) for c in coeffs[1:-1].split(",")])
        gen_name = special_names.get(label,'a')
        K = NumberField(poly, gen_name)
        nf_table[label] = K
    ffile.close()
    read_subfields()

def nf_lookup(label, verbose=False):
    r"""
    Returns a NumberField from its label, caching the result.
    """
    global nf_table, special_names
    if not nf_table:
        if verbose:
            print("reading field data from file")
        read_all_fields()
        if verbose:
            print(f"read {len(nf_table)} fields")
    if verbose:
        print(f"Looking up number field with label {label}")
    if label in nf_table:
        K = nf_table[label]
        if verbose: print(f"We have it: {K}")
        return K
    else:
        if verbose:
            print("We do not have it!")
        # for quadratic fields we can do this manually
        deg, s, absD, num = [ZZ(x) for x in label.split(".")]
        if deg==2:
            D = absD if s else -absD
            from sage.all import polygen
            x = polygen(QQ)
            pol = x**2 - x + (1-D)//4 if D%2 else x**2 - D//4
            K = NumberField(pol, 'a')
            assert K.discriminant()==D
            nf_table[label] = K
            if verbose:
                print(f"Adding quadratic field {K}")
            return K
        return None

def get_field_label(K, verbose=False, exact=True):
    r"""
    Return the label of field K (or None if it is not in nf_table)
    """
    global nf_table
    if not nf_table:
        if verbose:
            print("reading field data from file")
        read_all_fields()
        if verbose:
            print("read {} fields".format(len(nf_table)))
    if verbose:
        print("Looking up label of number field {}".format(K))
    if exact:
        return next((lab for lab in nf_table if K==nf_table[lab]), '')
    else:
        return next((lab for lab in nf_table if K.is_isomorphic(nf_table[lab])), '')

# from a field label return 'IQF', 'RQF', 'cubics', 'quartcs',
# quintics', or 'sextics':
def get_field_type_from_label(label):
    if label == '3.1.23.1':
        return 'gunnells'
    n,r,d,i = label.split(".")
    n = int(n) # degree
    r = int(r) # number of real embeddings
    if n==2:
        return ['IQF',None,'RQF'][r]
    return [None, None, None, 'cubics', 'quartics', 'quintics', 'sextics'][n]

# function to make a dict with field labels as keys and values lists
# of labels of subfields, excluding Q and the field itself.
def make_subfields():
    global subfield_table
    if subfield_table:
        return

    for label,K in nf_table.items():
        n = K.degree()
        subfield_table[label] = []
        if n not in [4,6]:
            continue
        for k,a,b in K.subfields():
            d = k.degree()
            if d==1 or d==n:
                continue
            sublabel = get_field_label(k, exact=False)
            assert sublabel # will fail if k not isomorphic to a field in the table
            subfield_table[label].append(sublabel)

def store_subfields(sffilename='ecnf_subfields'):
    make_subfields()
    with open(sffilename, 'w') as outfile:
        for label, sublabels in subfield_table.items():
            outfile.write(label)
            outfile.write(":")
            outfile.write(";".join(sublabels))
            outfile.write("\n")

def read_subfields(sffilename='ecnf_subfields'):
    global subfield_table
    subfield_table = {}
    with open(sffilename) as infile:
        for line in infile.readlines():
            k, subs = line.strip().split(":")
            subs = [] if subs=='' else subs.split(";")
            subfield_table[k] = subs

# Given a field label, return a list of the labels of its proper subfields (excluding Q)
def subfield_labels(label):
    if not subfield_table:
        read_all_fields() # also rads subfield data
    return subfield_table.get(label, [])


field_data = {} # dict whose keys are fields k and values are dicts holding:

# 'absD'                  # |disc(k)|
# 'G'                     # Gal(k/Q)
# 'autos'                 # automorphisms of k
# 'Plist'                 # sorted list of primes of k
# 'is_galois'             # is k Galois?
# 'label'                 # label of k
# 'poly'                  # defining polynomial of k as comma-sep. string
# 'nf_data'               # newform data for k (if available)
# 'used_curves'           # dict with key by norm(conductor), value a list of
#                              # curves so far processed witth that norm-conductor
# 'class_key'             # isogeny class sort key function (to be used with curves of the same
#                              # conductor only)
# 'cm_counts'             # dicts with keys conductor labels, values counts of
#                              # classes with rational CM (imaginary quadratic fields
#                              # only) for labelling of these, as they do not have
#                              # associated Bianchi newforms

# dict, keys are CM j-invariants, values are associated discriminanrs
cm_j_dict = {}

def add_field(K, field_label=None, prime_norm_bound=200, nf_data_file=None):
    global field_data, cm_j_dict
    if K in field_data:
        return

    if field_label==None:
        field_label = get_field_label(K)
        print("...created new label {}".format(field_label))
    print("Adding {} (label={}) to fields collection...".format(K,field_label))
    Kdata = {}
    Kdata['Plist'] = list(primes_iter(K,maxnorm=prime_norm_bound))
    absD = K.discriminant().abs()
    s = K.signature()[0] # number of real places
    d = K.degree()

# Warning: number fields whose label's 4'th component is not 1 will
# not be handled correctly here; this is only an issue when there is
# more than one field of given signature and abs(disc), so fine for
# quadratic fields.  This problem first hit for field 4.4.16448.2.
# When we processed the curves for that field they were in a different
# input file so were not mixed up with curves for 4.4.16448.1 luckily,
# and manual editing of the output was sufficient.

    Kdata['labels'] = field_label
    Kdata['absD'] = absD
    Kdata['autos'] = autos = K.automorphisms()
    Kdata['is_galois'] = len(autos)==K.degree()
    if d<5:
        Kdata['G'] = K.galois_group(names='b')
    for dd, f, j in cm_j_invariants_and_orders(K):
        cm_j_dict[j] = dd * (f ** 2)
    print("updated cm_j_dict") # " to {}".format(cm_j_dict))
    Kdata['used_curves'] = {}
    Kdata['nf_data'] = None
    Kdata['cm_counts'] = {}
    if nf_data_file:
        from files import read_newform_data
        assert d==2 and s==0 and absD in [3,4,7,8,11]
        print("reading newform data from {}".format(nf_data_file))
        Kdata['nf_data'] = read_newform_data(nf_data_file)
    field_data[K] = Kdata
    print("...finished adding field.")

import re
ideal_label_regex = re.compile(r'\d+\.\d+\.\d+')

def parse_ideal_label(s):
    if not ideal_label_regex.match(s):
        return [int(i) for i in s[1:-1].split(r',')]
    return [int(i) for i in s.split(r'.')]

def ideal_from_HNF(K,H):
    a,c,d = H
    return K.ideal([a,c+d*K.gen()])

def ideal_from_IQF_label(K,s):
    H = parse_ideal_label(s)
    H[0]/= H[2]
    return ideal_from_HNF(K,H)

def get_IQF_info(field_info_filename, maxpnorm=200, verbose=False):
    r"""
    Returns a number field and ordered list of primes.

    INPUT:

    - ``field_info_filename`` (string) -- name of data file.

    - ``maxpnorm`` (integer) -- bound on norm of primes to return.

    - ``verbose`` (boolean, default False) -- verbosity flag.

    OUTPUT:

    Tuple of a number field and a list of prime ideals, ordered as in the data file.
    """
    Plist=[]
    with open(field_info_filename) as field_info_file:
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
                P = ideal_from_IQF_label(K,lab)
                assert P.norm()==nm
                #print("Prime %s with char %s, degree %s, label %s, norm %s" % (P,p,deg,lab,nm))
                Plist.append(P)
    return K, Plist

