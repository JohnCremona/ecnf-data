import os
from sage.all import PolynomialRing, QQ
from files import ECNF_DIR
from fields import nf_lookup
from codec import NFelt

assert ECNF_DIR

def rewrite_ainvs(base_dir, field_label):
    """If we have created a curves file whose ainvs field looks like
    (a+1,a+1,a+1,-10*a-173,-27*a+899) instead of
    0,0;-1,0;1,0;-10,0;-20,0 then this function will read it and
    convert that field and rewrite it.
    """
    curves_filename_old = os.path.join(base_dir, 'curves.{}'.format(field_label))
    curves_filename_new = os.path.join(base_dir, 'curves.{}.new'.format(field_label))
    print("Reading from {}, writing to {}".format(curves_filename_old, curves_filename_new))
    R = PolynomialRing(QQ, 'a')
    K = nf_lookup(field_label)

    with open(curves_filename_old) as curves, open(curves_filename_new, 'w') as newcurves:
        for L in curves:
            data = L.split()
            ainvs = data[6]
            if ";" not in ainvs:
                print("Old line: {}".format(L))
                # convert to a list of 5 number field elts and hence to one string
                ainvs = ";".join([NFelt(K(R(a))) for a in ainvs[1:-1].split(",")])
                data[6] = ainvs
                L = " ".join(data)
                print("New line: {}\n".format(L))
            newcurves.write(L+"\n")
