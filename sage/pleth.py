# polynomial utilities using plethysms to compute "symmetric power" polynomials

from sage.all import QQ, prod, numerator, binomial, PolynomialRing, Sequence, SymmetricFunctions

def signed_coefficients(f):
    """Return list of coefficients of f with alternating signs, starting
    with signed leading coefficient and ending with constant term.

    """
    n = f.degree()
    return list(reversed([ci*(-1)**(n-i)
                          for i,ci in enumerate(f.coefficients(sparse=False))]))

def poly_from_signed_coefficients(c):
    """Return sum((-1)^i * c[i] * x^(n-i))
    """
    Rx = PolynomialRing(Sequence(c).universe(),'x')
    return Rx(list(reversed([ci*(-1)**i
                             for i,ci in enumerate(c)])))

def pleth_b(i,k,nmax=None):
    """Return the i,k plethysm, optionally limiting to partitions with no
    part greater than nmax, returning a list of (partition,
    multiplicity) pairs.
    """
    e = SymmetricFunctions(QQ).e()
    p = e[i].plethysm(e[k])
    if nmax:
        return [(part,mult)
                for part,mult in p if all(c<=nmax for c in part)]
    else:
        return list(p)

def ev(p, c):
    """Evaluate a list of (partition, multiplicity) pairs on a signed
    coefficients list c.
    """
    return sum([(prod(c[i] for i in part))*mult for part, mult in p])


def sym_pow(f,k):
    """Given a monic polynomial f of degree n, and k<=n, return the
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
        g = sym_pow(f,n-k)
        return numerator((-x)**binomial(n,k) * g(c/x) / c**binomial(n-1,k))

    N = binomial(n,k)
    cf = signed_coefficients(f)
    pp = [pleth_b(i,k,n) for i in range(N+1)]
    cc = [ev(p,cf) for p in pp]
    return poly_from_signed_coefficients(cc)
