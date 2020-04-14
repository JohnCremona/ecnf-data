# Billerey's algorithm for reducible primes, adapted from code by C. Schembri

from sage.all import polygen, ZZ, prod, primes, Set
from polys import star, star_iterate, r

def Frobenius_poly(E,qq):
    return E.reduction(qq).frobenius_polynomial()

def P_l(E, l):
    """Return Billerey's `P_l^*` as defined in (9).
    """
    P = polygen(ZZ)-1
    K = E.base_field()
    for q in K.primes_above(l):
        e = K(l).valuation(q)
        P = star(P,r(Frobenius_poly(E,q),12*e))
    return P

def B_l(E,l):
    """
    Return Billerey's `B_l` as defined after (9).
    """
    d = E.base_field().absolute_degree()
    P = P_l(E,l)
    return prod([ P(l**(12*k)) for k in range(1+d//2) ])

def R_q(E,q):
    """
    Return Billerey's R_q as defined in Theorem 2.8.
    """
    K = E.base_field()
    d = K.absolute_degree()
    h = K.class_number()
    P = r(Frobenius_poly(E,q),12*h)
    Q = r(((q**h).gens_reduced()[0]).absolute_minpoly(),12)
    return prod([ P.resultant(star_iterate(Q,k)) for k in range(1+d//2) ])

def B_bound(E, max_l=200, num_l=6, verbose=True):
    """Compute B_l for l up to max_l (at most) until num_l nonzero values
    are found (at most).  Return the list of primes dividing all B_l
    computed, excluding those dividing 6 or ramified or of bad
    reduction.  If no non-zero values are found return [0].
    """
    if verbose:
        print("Computing B-bound for {} with max_l={}, num_l={}".format(E.ainvs(),max_l,num_l))
    B = ZZ.zero()
    ells = []
    K = E.base_field()
    DK = K.discriminant()
    ED = E.discriminant().norm()
    B0 = ZZ(6*DK*ED)
    ll = primes(5,max_l) # iterator
    while len(ells)<num_l and B!=1:
        try:
            l = ll.next()
            while B0.valuation(l):
                l = ll.next()
        except StopIteration:
            break
        if verbose:
            print("..trying l={}".format(l))
        b = B_l(E,l)
        if b:
            if verbose:
                print("..ok, B_l = {}".format(b))
            if B:
                B = B.gcd(b)
            else:
                B = b.prime_to_m_part(B0)
            ells.append(l)
            if verbose:
                print("..so far, B = {} = {} using l in {}".format(B,B.support(),ells))

    if B:
        res = B.support()
        if verbose:
            print("..returning {}".format(res))
        return res
    # or we failed to find any nonzero values...
    if verbose:
        print("..failed to find a bound")
    return [0]

def R_bound(E, max_q=200, num_q=6, verbose=True):
    """Compute R_q for q dividing primes up to max_q (at most) until num_q
    nonzero values are found (at most).  Return the list of primes
    dividing all R_q computed, excluding those dividing 6 or ramified
    or of bad reduction.  If no non-zero values are found return [0].
    """
    if verbose:
        print("Computing R-bound for {} with max_q={}, num_q={}".format(E.ainvs(),max_q,num_q))
    B = ZZ.zero()
    ells = []
    K = E.base_field()
    DK = K.discriminant()
    ED = E.discriminant().norm()
    B0 = ZZ(6*DK*ED)
    ll = primes(5,max_q) # iterator
    while len(ells)<num_q and B!=1:
        try:
            l = ll.next()
            while B0.valuation(l):
                l = ll.next()
        except StopIteration:
            break
        q = K.prime_above(l)
        if verbose:
            print("..trying q={} above l={}".format(q,l))
        b = R_q(E,q)
        if b:
            if verbose:
                print("..ok, R_q = {}".format(b))
            if B:
                B = B.gcd(b)
            else:
                B = b.prime_to_m_part(B0)
            ells.append(l)
            if verbose:
                print("..so far, B = {} = {} using l in {}".format(B,B.support(),ells))

    if B:
        res = B.support()
        if verbose:
            print("..returning {}".format(res))
        return res
    # or we failed to find any nonzero values...
    if verbose:
        print("..failed to find a bound")
    return [0]

def Frobenius_filter(E,L, num_p=100):
    """Test whether E has an ell-isogeny locally.  Test if E mod P has an
    ell-isogeny for all P, using at most num_p primes P of good
    reduction.  If L is a list of primes, return the sublist which
    passes, otherwise L is a prime, and return True/False.

    Note that if (E, ell) fails the test then E certainly has no
    rational ell-isogeny, but not conversely since there are
    counterexamples to the local-global principal, but these are rare.
    """
    from sage.schemes.elliptic_curves.gal_reps_number_field import _maybe_borels
    return _maybe_borels(E,L,num_p)


def reducible_primes(E, num_l=6, max_l=200,
                        num_q=6, max_q=200, verbose=False):
    """Return a finite set of primes ell containing all those for which E
    has a K-rational ell-isogeny: i.e., the mod-ell representation is
    irreducible for all ell outside the set returned.

    We first compute Billeray's B_bound using at most num_l primes of
    size up to max_l.  If that fails we compute Billeray's R_bound
    using at most num_q primes of size up to max_q.  If that also
    fails we return [0].
    """
    if verbose:
        print("E = {}, seeking isogeny bound".format(E.ainvs()))
    K = E.base_field()
    DK = K.discriminant()
    ED = E.discriminant().norm()
    B0 = ZZ(6*DK*ED).prime_divisors()
    B1 = B_bound(E, max_l, num_l, verbose)
    if B1 == [0]:
        print("...  B_bound ineffective using max_l={}, moving on to R-bound".format(max_l))
        B1 = R_bound(E, max_q, num_q, verbose)
        if B1 == [0]:
            print("... R_bound ineffective using max_q={}",format(max_q))
            return [0]
        if verbose:
            print("... R_bound = {}".format(B1))
    else:
        if verbose:
            print("... B_bound = {}".format(B1))
    B = sorted(Set(B0 + B1))
    if verbose:
        print("... combined bound = {}".format(B))
    B = Frobenius_filter(E,B)
    if verbose:
        print("... after Frobenius filter = {}".format(B))
    return B
