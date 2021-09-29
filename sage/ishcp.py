from sage.all import GF, polygen, ceil, floor, ZZ, pari, prime_range, PolynomialRing
from sage.libs.pari.convert_sage import gen_to_sage
import six
if six.PY3:
    from functools import reduce
# in python2 reduce was builtin but not in python3

def SupersingularPolynomial(p):
    """ Implementation of Finotti's forumla for the supersingular polynomial """
    F = GF(p)
    def prod_p(list):
        return reduce(lambda x,y: F(x)*F(y), list, F(1))
    def binomial_p(r,s):
        assert s >= 0
        if r < s:
            return F(0)
        if s < r-s:
            s = r-s
        return prod_p(range(s+1,r+1)) / prod_p(range(1,r-s+1))
    x = polygen(F)
    if p <= 5:
        return x
    r = (p-1) // 2
    r1, r2 = ceil(r/3), floor(r/2)
    s1, s2 = floor(r/3), ceil(r/2)
    c, c1728 = F(-27)/F(4), F(1728)
    f = F(-2)**r*sum([binomial_p(r,i)*binomial_p(i,3*i-r)*c**i*x**(i-s1)*(x-c1728)**(s2-i) for i in range(r1,r2+1)])
    return f


def IsHilbertClassPolynomial(f):
    """ Returns true if f is a Hilbert class polynomial, false otherwise. """
    T = polygen(ZZ)
    X,Y = PolynomialRing(ZZ,2).gens()
    f = f.change_ring(ZZ)(T)
    assert f.is_monic() and f.is_irreducible()
    Fq = GF(101)
    Tp = polygen(Fq)
    Xp,Yp = PolynomialRing(Fq,2).gens()
    for p in prime_range(1,100):
        phi = gen_to_sage(pari.polmodular(p),{'x':X,'y':Y})
        # perform resultant check in Fq first (we expect it to fail)
        if f(Tp).divides(phi([Xp,Yp]).resultant(f(Yp),Yp)([Tp,0])):
            res = phi.resultant(f(Y),Y)([T,0])
            if f.divides(res):
                if (f**2).divides(res):
                    return True, p
                else:
                    continue
        if p <= 3:
            continue
        t = polygen(GF(p))
        g = f(t)
        if g.gcd(g.derivative()) != 1:
            continue
        if not g.divides(SupersingularPolynomial(p)(t)):
            return False, p
    assert False, "p > 100, this really should never happen"
