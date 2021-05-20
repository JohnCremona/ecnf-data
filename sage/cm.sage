from sage.all import QuadraticField, NumberField, srange, EllipticCurve, polygen, QQ, is_squarefree
from codec import ideal_to_string
from psort import ideal_label, primes_iter
from nfscripts import ap_list #make_curves_line, make_ec_dict

x = polygen(QQ)

def IQF(d):
    t,n = (0,d) if d%4 in [1,2] else (1,(d+1)//4)
    return NumberField(x^2-t*x+n, 'a')


Dfj_list = cm_j_invariants_and_orders(QQ)

cm_j_invs = dict([(D.abs() if D%2 else -D//4, dict([(f,j) for d,f,j in Dfj_list if d==D]))
                  for D in set(D for D,f,j in Dfj_list)])

d_list = sorted([d for d in cm_j_invs])
fields = dict([(d,IQF(d)) for d in d_list])


def K_even_iterator(K, norm_bound, S=[]):
    n = K.discriminant().abs()//4
    By = (norm_bound//n).isqrt()
    for y in srange(By+1):
        Bx = (norm_bound-n*y*y).isqrt()
        xmin = 1
        if y: xmin = -Bx
        for x in srange(xmin,Bx+1):
            a = K([x,y])
            if all([a.valuation(p)==0 for p in S]):
                yield a

def K_odd_iterator(K, norm_bound, S=[]): # for all except K1, K2
    n = K.discriminant().abs()
    By = (4*norm_bound//n).isqrt()
    for y in srange(By+1):
        Bz = (4*norm_bound-n*y*y).isqrt()
        zmin = y+2
        for z in srange(zmin,Bz+1,2):
            x = (z-y)/2
            a = K([x,y])
            if all([a.valuation(p)==0 for p in S]):
                yield a

def K_iterator(K, norm_bound, S=[]):
    D = K.discriminant().abs()
    if D%4==0:
        return K_even_iterator(K, norm_bound, S)
    else:
        return K_odd_iterator(K, norm_bound, S)

def is_powerfree(a, d=2):
    return all([e<d for p,e in a.factor()])

def all_non_iso(EE):
    return all([not any([EE[j].is_isomorphic(EE[k]) for k in range(j+1,len(EE))]) for j in range(len(EE))])

def uniq_iso(EE):
    return [E for i,E in enumerate(EE) if not any([E.is_isomorphic(EE[j]) for j in range(i)])]


def curves_K1(max_norm, min_norm=1, f=1, verb=False):
    K = fields[1]
    j = cm_j_invs[1][f]
    m = [0,4,2][f]
    E = EllipticCurve(j=j).change_ring(K)
    S = K(6).support()
    supps = {}
    twists = {}
    for d1 in K_iterator(K, max_norm.isqrt(),S):
        if is_powerfree(d1,2):
            print("d1={}".format(d1))
            d1supp = d1.support()
            t = tuple(d1supp)
            if t in supps:
                continue
            supps[t] = 1
            for d in K.selmer_group_iterator(S+d1supp,m):
                if d in twists:
                    continue
                twists[d] = 1
                #print("d={}".format(d))
                if m==2:
                    Ed = E.quadratic_twist(d)
                else:
                    Ed = E.quartic_twist(d)
                Nn = Ed.conductor().norm()
                #print("cond.norm={}".format(Nn))
                if min_norm <= Nn <= max_norm:
                    Ed = Ed.global_minimal_model()
                    if verb: print(d, Ed.ainvs(), Nn)
                    yield Ed

def curves_K2(max_norm, min_norm=1, verb=False):
    K = fields[2]
    j = cm_j_invs[2][1] # = 8000
    E = EllipticCurve(j=j).change_ring(K)
    S = K(6).support()
    for d1 in K_iterator(K, max_norm.isqrt(),S):
        if is_squarefree(d1):
            for d0 in K.selmer_group_iterator(S,2):
                d = d0*d1
                Ed = E.quadratic_twist(d)
                Nn = Ed.conductor().norm()
                if min_norm <= Nn <= max_norm:
                    Ed = Ed.global_minimal_model()
                    if verb: print(d, Ed.ainvs(), Nn)
                    yield Ed


def curves_K3(max_norm, min_norm=1, f=1, verb=False):
    K = fields[3]
    j = cm_j_invs[3][f]
    m = [0,6,2,2][f]
    E = EllipticCurve(j=j).change_ring(K)
    S = K(6).support()
    supps = {}
    twists = {}
    for d1 in K_iterator(max_norm.isqrt(),S):
        if is_powerfree(d1,2):
            print("d1={}".format(d1))
            d1supp = d1.support()
            t = tuple(d1supp)
            if t in supps:
                continue
            supps[t] = 1
            for d in K.selmer_group_iterator(S+d1supp,m):
                if d in twists:
                    continue
                twists[d] = 1
                #print("d={}".format(d))
                if m==2:
                    Ed = E.quadratic_twist(d)
                else:
                    Ed = E.sextic_twist(d)
                Nn = Ed.conductor().norm()
                if min_norm <= Nn <= max_norm:
                    Ed = Ed.global_minimal_model()
                    if verb: print(d, Ed.ainvs(), Nn)
                    yield Ed

def curves_K7(max_norm, min_norm=1, f=1, verb=False):
    K = fields[7]
    j = cm_j_invs[7][f]
    E = EllipticCurve(j=j).change_ring(K)
    S = K(6*7).support()
    for d1 in K_iterator(K, max_norm.isqrt(),S):
        if is_squarefree(d1):
            for d0 in K.selmer_group_iterator(S,2):
                d = d0*d1
                Ed = E.quadratic_twist(d)
                Nn = Ed.conductor().norm()
                if min_norm <= Nn <= max_norm:
                    Ed = Ed.global_minimal_model()
                    if verb: print(d, Ed.ainvs(), Nn)
                    yield Ed

# use for K from 11 on (odd, only one f, only 2 units)
def curves_K(d, max_norm, min_norm=1, verb=False):
    K = fields[d]
    j = cm_j_invs[d][1]
    E = EllipticCurve(j=j).change_ring(K)
    S = K(6*11).support()
    for d1 in K_iterator(K, max_norm.isqrt(),S):
        if is_squarefree(d1):
            for d0 in K.selmer_group_iterator(S,2):
                d = d0*d1
                Ed = E.quadratic_twist(d)
                Nn = Ed.conductor().norm()
                if min_norm <= Nn <= max_norm:
                    Ed = Ed.global_minimal_model()
                    if verb: print(d, Ed.ainvs(), Nn)
                    yield Ed


def dump(Elist, outfilename=None):
    if outfilename == None:
        return
    outfile = open(outfilename, mode='w')
    K = Elist[0].base_field()
    Plist = list(primes_iter(K,maxnorm=100))
    def sort_key(E):
        return (E.conductor().norm(), E.conductor(), ap_list(E,Plist))

    field_lab = "2.0.{}.1".format(K.discriminant().abs())
    outfile.write("Field {}\n".format(field_lab))
    cond_lab_count = {}
    print("Before sorting:")
    for E in Elist:
        print(E.ainvs())
    Elist.sort(key=sort_key)
    print("After sorting:")
    for E in Elist:
        print(E.conductor().norm(), E.conductor(), E.ainvs(), ap_list(E, Plist[:5]))

    for E in Elist:
        cond = E.conductor()
        cond_label = ideal_label(cond)
        class_label_base = "{}-CM".format(cond_label)
        i = cond_lab_count.get(class_label_base,0) + 1
        cond_lab_count[class_label_base] = i
        class_label = class_label_base + 'abcdef'[i-1]
        outfile.write("Conductor {}\n".format(ideal_to_string(cond)))
        outfile.write("Isogeny_class {}\n".format(class_label))
        ainvs = str(list(E.ainvs())).replace("a", "w")
        outfile.write("Curve {}\n".format(ainvs))
    outfile.close()


def cm_curves(field,max_norm, min_norm=1, outfilename=None, verbose=False):
    """ field is 1, 2, 3, 7, 11, 19, 43, 67, 163
    """
    if field==1:
        if verbose:
            print("Qsqrt-1, norms {}-{}".format(min_norm,max_norm))
        E1f1 = list(curves_K1(max_norm,min_norm, 1, verbose))
        E1f2 = list(curves_K1(max_norm,min_norm, 2, verbose))
        E1 = E1f1 + E1f2
        assert all_non_iso(E1)
        if verbose:
            print(" found {}+{}={} curves".format(len(E1f1),len(E1f2),len(E1)))
        dump(E1, outfilename)
        return E1

    if field==2:
        if verbose:
            print("Qsqrt-2, norms {}-{}".format(min_norm,max_norm))
        E2 = list(curves_K2(max_norm, min_norm,verbose))
        assert all_non_iso(E2)
        if verbose:
            print(" found {} curves".format(len(E2)))
        dump(E2, outfilename)
        return E2

    if field==3:
        if verbose:
            print("Qsqrt-3, norms {}-{}".format(min_norm,max_norm))
        E3f1 = list(curves_K3(max_norm, min_norm,1,verbose))
        E3f2 = list(curves_K3(max_norm, min_norm,2,verbose))
        E3f3 = list(curves_K3(max_norm, min_norm,3,verbose))
        E3 = E3f1+E3f2+E3f3
        assert all_non_iso(E3)
        if verbose:
            print(" found {}+{}+{}={} curves".format(len(E3f1),len(E3f2),len(E3f3),len(E3)))
        dump(E3, outfilename)
        return E3

    if field==7:
        if verbose:
            print("Qsqrt-7, norms {}-{}".format(min_norm,max_norm))
        E7f1 = list(curves_K7(max_norm, min_norm,1,verbose))
        assert all_non_iso(E7f1)
        E7f2 = list(curves_K7(max_norm, min_norm,2,verbose))
        assert all_non_iso(E7f2)
        E7 = E7f1+E7f2
        assert all_non_iso(E7)
        if verbose:
            print(" found {}+{}={} curves".format(len(E7f1),len(E7f2),len(E7)))
        dump(E7, outfilename)
        return E7

    # now field >=11
    if verbose:
        print("Qsqrt-{}, norms {}-{}".format(field, min_norm,max_norm))

    EE = list(curves_K(field, max_norm, min_norm,verbose))
    assert all_non_iso(EE)
    if verbose:
        print(" found {} curves".format(len(EE)))
    dump(EE, outfilename)
    for E in EE:
        print(E.ainvs(), E.conductor(), E.conductor().norm())
    return EE

