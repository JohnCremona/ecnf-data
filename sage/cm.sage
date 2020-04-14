from sage.all import QuadraticField, NumberField, srange, EllipticCurve, polygen, QQ, is_squarefree
from nfscripts import make_curves_line, make_ec_dict

x = polygen(QQ)
K1 = QuadraticField(-1, 'i'); i=K1.gen()
K2 = QuadraticField(-2, 't'); t=K2.gen()
K3 = NumberField(x^2-x+1, 'z'); z=K3.gen()
K7 = NumberField(x^2-x+2, 'a'); a=K7.gen()
K11 = NumberField(x^2-x+3, 'b'); b=K11.gen()

def K1_iterator(norm_bound, S=[]):
    By = norm_bound.isqrt()
    for y in srange(By+1):
        Bx = (norm_bound-y*y).isqrt()
        xmin = 1
        for x in srange(xmin,Bx+1):
            a = x+y*i
            if all([a.valuation(p)==0 for p in S]):
                yield a

def K2_iterator(norm_bound, S=[]):
    By = (norm_bound//2).isqrt()
    for y in srange(By+1):
        Bx = (norm_bound-2*y*y).isqrt()
        xmin = 1
        if y: xmin = -Bx
        for x in srange(xmin,Bx+1):
            a = x+y*t
            if all([a.valuation(p)==0 for p in S]):
                yield a

def K3_iterator(norm_bound, S=[]):
    By = (4*norm_bound//3).isqrt()
    for y in srange(By+1):
        Bz = (4*norm_bound-3*y*y).isqrt()
        zmin = y+2
        for z in srange(zmin,Bz+1,2):
            x = (z-y)/2
            a = x+y*K3.gen()
            if all([a.valuation(p)==0 for p in S]):
                yield a

def K7_iterator(norm_bound, S=[]):
    By = (4*norm_bound//7).isqrt()
    for y in srange(By+1):
        Bz = (4*norm_bound-7*y*y).isqrt()
        zmin = y+2
        if y:
            zmin = -Bz
            if (y-zmin)%2: zmin +=1
        for z in srange(zmin,Bz+1,2):
            x = (z-y)/2
            a = x+y*K7.gen()
            if all([a.valuation(p)==0 for p in S]):
                yield a

def K11_iterator(norm_bound, S=[]):
    By = (4*norm_bound//11).isqrt()
    for y in srange(By+1):
        Bz = (4*norm_bound-11*y*y).isqrt()
        zmin = y+2
        if y:
            zmin = -Bz
            if (y-zmin)%2: zmin +=1
        for z in srange(zmin,Bz+1,2):
            x = (z-y)/2
            a = x+y*K11.gen()
            if all([a.valuation(p)==0 for p in S]):
                yield a

def is_powerfree(a, d=2):
    return all([e<d for p,e in a.factor()])

def all_non_iso(EE):
    return all([not any([EE[j].is_isomorphic(EE[k]) for k in range(j+1,len(EE))]) for j in range(len(EE))])

def uniq_iso(EE):
    return [E for i,E in enumerate(EE) if not any([E.is_isomorphic(EE[j]) for j in range(i)])]



def curves_K1(max_norm, min_norm=1, f=1, verb=False):
    j = [0,1728,287496][f]
    m = [0,4,2][f]
    E = EllipticCurve(j=j).change_ring(K1)
    S = K1(6).support()
    supps = {}
    twists = {}
    for d1 in K1_iterator(max_norm.isqrt(),S):
        if is_powerfree(d1,2):
            print("d1={}".format(d1))
            d1supp = d1.support()
            t = tuple(d1supp)
            if t in supps:
                continue
            supps[t] = 1
            for d in K1.selmer_group_iterator(S+d1supp,m):
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
    j = 8000
    E = EllipticCurve(j=j).change_ring(K2)
    S = K2(6).support()
    for d1 in K2_iterator(max_norm.isqrt(),S):
        if is_squarefree(d1):
            for d0 in K2.selmer_group_iterator(S,2):
                d = d0*d1
                Ed = E.quadratic_twist(d)
                Nn = Ed.conductor().norm()
                if min_norm <= Nn <= max_norm:
                    Ed = Ed.global_minimal_model()
                    if verb: print(d, Ed.ainvs(), Nn)
                    yield Ed

def curves_K7(max_norm, min_norm=1, f=1, verb=False):
    j = [0,-3375,16581375][f]
    E = EllipticCurve(j=j).change_ring(K7)
    S = K7(6*7).support()
    for d1 in K7_iterator(max_norm.isqrt(),S):
        if is_squarefree(d1):
            for d0 in K7.selmer_group_iterator(S,2):
                d = d0*d1
                Ed = E.quadratic_twist(d)
                Nn = Ed.conductor().norm()
                if min_norm <= Nn <= max_norm:
                    Ed = Ed.global_minimal_model()
                    if verb: print(d, Ed.ainvs(), Nn)
                    yield Ed


def curves_K11(max_norm, min_norm=1, verb=False):
    j = -32768
    E = EllipticCurve(j=j).change_ring(K11)
    S = K11(6*11).support()
    for d1 in K11_iterator(max_norm.isqrt(),S):
        if is_squarefree(d1):
            for d0 in K11.selmer_group_iterator(S,2):
                d = d0*d1
                Ed = E.quadratic_twist(d)
                Nn = Ed.conductor().norm()
                if min_norm <= Nn <= max_norm:
                    Ed = Ed.global_minimal_model()
                    if verb: print(d, Ed.ainvs(), Nn)
                    yield Ed


def curves_K3(max_norm, min_norm=1, f=1, verb=False):
    j = [0,0,54000,-12288000][f]
    m = [0,6,2,2][f]
    E = EllipticCurve(j=j).change_ring(K3)
    S = K3(6).support()
    supps = {}
    twists = {}
    for d1 in K3_iterator(max_norm.isqrt(),S):
        if is_powerfree(d1,2):
            print("d1={}".format(d1))
            d1supp = d1.support()
            t = tuple(d1supp)
            if t in supps:
                continue
            supps[t] = 1
            for d in K3.selmer_group_iterator(S+d1supp,m):
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

def dump(Elist, outfilename=None):
    if outfilename == None:
        return
    outfile = open(outfilename, mode='w')
    for E in Elist:
        outfile.write(make_curves_line(make_ec_dict(E))+'\n')
    outfile.close()


def cm_curves(field,max_norm, min_norm=1, outfilename=None, verbose=False):
    """ field is 1, 2, 3, 7, 11
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

    if field==11:
        if verbose:
            print("Qsqrt-11, norms {}-{}".format(min_norm,max_norm))

    E11 = list(curves_K11(max_norm, min_norm,verbose))
    assert all_non_iso(E11)
    if verbose:
        print(" found {} curves".format(len(E11)))
    dump(E11, outfilename)
    return E11

