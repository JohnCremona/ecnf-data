from sage.rings.all import operator

def star(P,Q):
    """P*Q where for P,Q monic polynomials with distinct roots a_i, b_j
    we set P*Q to be the monic polynomial of degree deg(P)*deg(Q) with
    roots a_i*b_j
    """
    return P.composed_op(Q, operator.mul, monic=True)

def r(P, n):
    """Return the polynomial whose roots are the n'th powers of those of P
    """
    return P.power_roots(n,monic=True)

def star_iterate(P, k):
    """Return k-fold star-iterate P*P*...*P, for k>=0.

    We ought to be able to use sage.groups.generic.multiple but cannot
    since for a general op that requires one to provide an inverse op,
    even if k>=0, which is silly.
    """
    return P.compose_power(k, monic=True)

# The following is used in my efficiency improvement to Larsen's algorithm

def s(f,k):
    """Given f monic of degree n with distinct roots, returns the monic
    polynomial of degree binomial(n,k) whose roots are all products of
    k distinct roots of f.
    """
    n = f.degree()
    x = f.variables()[0]
    if k==0:
        return x-1
    if k==1:
        return f
    c = (-1)**n * f(0)
    if k==n:
        return x-c
    if k>n-k: # use s(f,n-k)
        g = s(f,n-k)
        from sage.arith.all import binomial
        return ((-x)**binomial(n,k) * g(c/x) / c**binomial(n-1,k)).numerator()

    if k==2:
        return (star(f,f)//r(f,2)).nth_root(2)
    if k==3:
        f2 = s(f,2)
        return (star(f2,f)*r(f,3) // star(r(f,2),f)).nth_root(3)

    fkn = fkd = 1
    for j in range(1,k+1):
        g = star(r(f,j), s(f,k-j))
        if j%2:
            fkn *= g
        else:
            fkd *= g

    fk = fkn//fkd
    assert fk*fkd==fkn
    return fk.nth_root(k)

