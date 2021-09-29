from sage.all import ZZ, Primes, FiniteField, EllipticCurve, gcd, hilbert_class_polynomial

def cm_test(f, maxp=1000):
    if not f.is_monic():
        return False
    if not all([c in ZZ for c in f.coefficients()]):
        return False
    bads = [6*f.discriminant(), f(0), f(1728)]
    cmd = 0
    cmf = 0
    for p in Primes():
        if p>maxp:
            break
        if any([p.divides(b) for b in bads]):
            continue
        glist = [g for g,e in f.factor_mod(p)]
        assert all([g.is_irreducible() for g in glist])
        # Flist = [FiniteField(p**g.degree(), name='j', modulus=g) for g in glist]
        # Elist = [EllipticCurve(F, j=F.gen()) for F in Flist]
        Flist = [FiniteField(p**g.degree()) for g in glist]
        Elist = [EllipticCurve(F, j=g.roots(F)[0][0]) for g,F in zip(glist,Flist)]
        tlist = [E.trace_of_frobenius() for E in Elist]
        ss = [t%p for t in tlist]
        consistent = all(ss) or not any(ss)
        if not consistent:
            return False
        if not ss[0]:
            continue
        dlist = [t**2-4*p**g.degree() for g,t in zip(glist, tlist)]
        d0list = [d.squarefree_part() for d in dlist]
        if not all([d0==d0list[0] for d0 in d0list]):
            return False
        D0 = d0list[0]
        ff = gcd([(d//D0).isqrt() for d in dlist])
        if cmd==0:
            cmd = D0
            cmf = ff
        elif cmd != D0:
            return False
        else:
            cmf = cmf.gcd(ff)
        #print("so far D0={}, f={}".format(cmd,cmf))

    if cmd % 4 != 1:
        cmd = cmd * 4
        assert cmf%2==0
        cmf = cmf // 2

    h = f.degree()
    for ff in cmf.divisors():
        d = cmd*ff**2
        if d.class_number() == h:
            if hilbert_class_polynomial(d) == f:
                return True, (cmd, ff)
    return False

