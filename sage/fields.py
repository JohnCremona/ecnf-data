# We have to look up number fields in the database from their labels,
# but only want to do this once for each label, so we will maintain a
# dict of label:field pairs:

from sage.all import QQ, ZZ, PolynomialRing, NumberField
nf_table = {}

special_names = {'2.0.4.1': 'i',
                 '2.2.5.1': 'phi',
                 '4.0.125.1': 'zeta5',
                 }

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
            print("read {} fields".format(len(nf_table)))
    if verbose:
        print("Looking up number field with label {}".format(label))
    if label in nf_table:
        K = nf_table[label]
        if verbose: print("We have it: {}".format(K))
        return K
    else:
        if verbose:
            print("We do not have it!")
        return None
