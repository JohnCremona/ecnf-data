# Check that the CM discriminants in curves files are correct

import os
from sage.all import cm_j_invariants_and_orders
from fields import nf_lookup
from codec import parse_curves_line, encoders
from files import ECNF_DIR, all_ftypes

cm_data = {}

def add_field_cm_data(field_label):
    global cm_data
    if field_label in cm_data:
        return
    K = nf_lookup(field_label)
    def fix(D):
        return abs(D) if K(D).is_square() else D
    cm_data[field_label] = dict([(encoders['jinv'](j), fix(d*f**2)) for d, f, j in cm_j_invariants_and_orders(K)])
    print("updated cm_data for {} with {}".format(field_label, cm_data[field_label]))
    return

def check_one_cm(field_label, j_string, CMD, verbose=False):
    add_field_cm_data(field_label)
    if verbose:
        print("Checking CM discriminant {} for {}, {}:".format(CMD, field_label, j_string))
    D = cm_data[field_label].get(j_string, 0)
    if D==CMD:
        return True, D
    else:
        print("Field {}, j = {}, discriminant was [{}] but should be [{}]".format(field_label, j_string, CMD, D))
        return False, D

def check_field_cm(ftype, field_label, verbose=False):
    cfile = os.path.join(ECNF_DIR, ftype, "curves.{}".format(field_label))
    print("Reading from {}".format(cfile))
    n = n_good = n_bad = 0
    with open(cfile) as curves:
        for L in curves:
            label, record = parse_curves_line(L)
            if label:
                n += 1
                ok, D = check_one_cm(field_label, record['jinv'], record['cm'], verbose)
                if ok:
                    n_good +=1
                else:
                    n_bad +=1
    print("Read {} curves from {}: {} are good, {} are bad".format(n,cfile, n_good, n_bad))
    return (n_bad==0)

def check_cm(ftypes=all_ftypes, verbose=False):
    bad_fields = []
    for ft in ftypes:
        if ft=="sengun":
            print("ignoring curves in {}".format(ft))
            continue
        print("checking curves in {}".format(ft))
        with open(os.path.join(ECNF_DIR,ft,"fields.txt")) as field_file:
            flds = [f[:-1] for f in field_file.readlines()]
        for f in flds:
            print("checking curves over {}".format(f))
            ok = check_field_cm(ft, f, verbose=verbose)
            if not ok:
                bad_fields.append(f)
    if bad_fields:
        print("Fields with bad cm columns: {}".format(bad_fields))
    else:
        print("All cm columns are correct")
