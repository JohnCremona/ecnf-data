# Code to check that a curve model is global minimal or semi-global minimal

import sys
from files import read_curves_new

def check_curve_model(E):
    """
    E is a curve defined over a number field
    """
    Emin = E.global_minimal_model()  if E.base_field().class_number()==1 else E.global_minimal_model(semi_global=True)
    return E.discriminant().norm() == Emin.discriminant().norm()

def check_curves(infile, verbose=False):
    bad_curves = {}
    n = 0
    nbad = 0
    for (field_label,N_label,N_def,iso_label,c_num,E) in read_curves_new(infile):
        n += 1
        if not check_curve_model(E):
            nbad += 1
            label = "{}-{}-{}{}".format(field_label, N_label, iso_label, c_num)
            print("{} is bad".format(label))
            bad_curves[label] = E
        if n%100==0 and verbose:
            print("{}...".format(n), end="")
            sys.stdout.flush()
            if n%1000==0:
                print()
    if nbads:
        ok = "!!"
    else:
        ok = "OK"
    print("{}: {} bad curves out of {}".format(ok,nbad,n))
    return bad_curves

def check_type(t):
    fields = read_data("../{}/fields.txt".format(t), str)
    all_bads = {}
    for f in fields:
        print(f)
        all_bads[f] = check_curves("../{}/curves.{}".format(t,f))

    return all_bads
