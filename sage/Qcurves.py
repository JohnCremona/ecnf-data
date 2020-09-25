# Qcurves.py: Implementation of a test for elliptic curves over
# arbitrary number fields to be Q-curves, as described in Cremona and
# Najman, "Q-curves over odd degree number fields",
# https://arxiv.org/abs/2004.10054

#######################################################################
#
# Copyright 2020 John Cremona
#
# This is free software; you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the
# Free Software Foundation; either version 2 of the License, or (at your
# option) any later version.
#
# This code is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
# FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
# for more details.
#
# You should have received a copy of the GNU General Public License
# along with this file; if not, write to the Free Software Foundation,
# Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301, USA
#
#######################################################################


from sage.all import (ZZ, QQ, EllipticCurve, primes, hilbert_class_polynomial, cm_j_invariants, NumberField, pari, polygen, LCM, Set )
from sage.schemes.elliptic_curves.cm import discriminants_with_bounded_class_number

# To avoid recomputing the CM j-invariants of any fixed degree, we
# keep a dict of all the j-polynomials in each degree.

# The following lines initialise the dict unless it has already
# started to be filled (so we can reload this file without losing
# values already computed):

try:
    assert len(j_pols_by_degree)
except:
    j_pols_by_degree = {}

def is_cm_pol(H):
    """H is a minimal polynomial of an algebraic number, hence
    irreducible and monic but not necessarily integral.

    Returns either a negative discriminant D, if H is the minimal
    polynomial of a CM j-invariant with discriminant D, or 0
    otherwise.

    Uses a cache, so that for each degree d we compute the full set of
    CM polynomials of degree d only once.

    """
    global j_pols_by_degree
    d = H.degree()
    if d in j_pols_by_degree:
        if H in j_pols_by_degree[d]:
            return j_pols_by_degree[d][H]
        return ZZ(0)

    if not H.is_monic():
        return ZZ(0)
    if not all(c in ZZ for c in H.coefficients()):
        return ZZ(0)

    print("Computing all degree {} CM j-polynomials".format(d))
    Dlist = [D*f*f for D,f in discriminants_with_bounded_class_number(d)[d]]
    j_pols_by_degree[d] = dict([(hilbert_class_polynomial(D),D) for D in Dlist])
    print(" done: {} discriminants {}".format(len(Dlist), Dlist))

    if H in j_pols_by_degree[d]:
        return j_pols_by_degree[d][H]
    return ZZ(0)

def is_Q_curve(E, maxp=100, certificate=False, verbose=False):
    r"""Return True if E is a Q-curve, i.e. is isogenous (over Qbar) to all
    its Galois conjugates.

    ``maxp`` (int, default 100): bound on primes used for checking
    necessary local conditions.  The result will not depend on this,
    but using a larger value is more likely to return False faster.

    ``certificate`` (bool, default ``False``): if ``True`` then a
    second value is returned giving a certificate for the Q-curve
    property.  This a dict, either {'CM':D} where D is a negative
    discriminant, or {'CM':0, 'core_poly', 'rho', 'r', 'N'} where
    'core_poly' is an irreducible monic polynomial over `QQ` of degree
    `2**\rho`, all of whose roots are j-invariants of curves isogenous
    to $E$, `N` (the core level) is a square-free integer with `r`
    prime factors which is the LCM of the degrees of the isogenies
    between these conjugates.  For example, if `j(E')=j\in\QQ` for any
    curve $E'$ isogenous to $E$, then the certificate is {'CM':0,
    'r':0, 'rho':0, 'core_poly': x-j, 'N':1}.

    Method:

    1. If E has rational j-invariant or has CM then True.

    2. Replace E by a curve defined over K=Q(j). Let N be the conductor norm.

    3. For all p|N check that the valuations of j at all P|p are
    either all negative or all non-negative; else False.

    4. For p<=maxp not dividing N, check that either all P|p are
    ordinary or all are supersingular, else False.  If all are
    ordinary, check that the integers a_P(E)^2-4*N(P) have the same
    square-free part; else False.

    5. Compute the K-isogeny class of E using the "heuristic" option
    (so it is not guaranteed to be complete).  Check whether the set
    of j-invariants of curves in the class of 2-power degree contains
    a complete Galois orbit.  If so:True.

    6. Otherwise repeat step 4 for more primes, and if still
    undecided, repeat Step 5 without the "heuristic" option, to get
    the complete K-isogeny class (which will probably be no bigger
    than before).  Now return True if the set of j-invariants of
    curves in the class contains a complete Galois orbit, otherwise
    return False.

    """
    if verbose:
        print("Checking whether {} is a Q-curve".format(E))

    # Step 1

    # all curves with rational j-invariant are Q-curves:
    jE = E.j_invariant()
    if jE in QQ:
        if verbose:
            print("Yes: j(E) is in QQ")
        if certificate:
            x = polygen(QQ)
            if jE in cm_j_invariants(QQ):
                return True, {'CM': is_cm_pol(x-jE)}
            else:
                return True, {'CM': ZZ(0), 'r': ZZ(0), 'rho': ZZ(0), 'N': ZZ(1), 'core_poly': x}
        else:
            return True

    # CM curves are Q-curves:
    jpoly = jE.minpoly()
    D = is_cm_pol(jpoly)
    if D:
        if verbose:
            print("Yes: E is CM (discriminant {})".format(D))
        if certificate:
            return True, {'CM': D}
        else:
            return True

    # Step 2: replace E by a curve defined over Q(j(E)):

    K = E.base_field()
    if jpoly.degree()<K.degree():
        if verbose:
            print("switching to smaller base field: j's minpoly is {}".format(jpoly))
        f = pari(jpoly).polredbest().sage({'x':jpoly.parent().gen()})
        K2 = NumberField(f, 'b')
        jE = jpoly.roots(K2)[0][0]
        if verbose:
            print("New j is {} over {}, with minpoly {}".format(jE, K2, jE.minpoly()))
        assert jE.minpoly()==jpoly
        E = EllipticCurve(j=jE)
        K = K2
        if verbose:
            print("New test curve is {}".format(E))

    # Step 3: check primes of bad reduction

    NN = E.conductor().norm()
    for p in NN.support():
        Plist = K.primes_above(p)
        if len(Plist)<2:
            continue
        pot_mult = [jE.valuation(P)<0 for P in Plist]
        consistent = all(pot_mult) or not any(pot_mult)
        if not consistent:
            if verbose:
                print("No: inconsistency at the {} primes dividing {} ".format(len(Plist),p))
                print("  - potentially multiplicative: {}".format(pot_mult))
            if certificate:
                return False, None
            else:
                return False

    # Step 4 check: primes P of good reduction above p<=B:

    if verbose:
        print("Applying local tests at good primes above p<={}".format(maxp))

    res4, p = Step4Test(E, B=maxp, oldB=0, verbose=verbose)
    if not res4:
        if verbose:
            print("No: local test at p={} failed".format(p))
        if certificate:
            return False, None
        else:
            return False

    if verbose:
        print("...all local tests pass for p<={}".format(maxp))

    # Step 5: compute the (partial) K-isogeny class of E and test the
    # set of j-invariants in the class:

    C = E.isogeny_class(algorithm='heuristic', minimal_models=False)
    jC = [E2.j_invariant() for E2 in C]
    res, centrejpols = conjugacy_test(jC, verbose=verbose)
    if res:
        if verbose:
            print("Yes: the isogeny class contains a complete conjugacy class of j-invariants")
        if certificate:
            for f in centrejpols:
                rho = f.degree().valuation(2)
                centre_indices = [i for i,j in enumerate(jC) if f(j)==0]
                M = C.matrix()
                core_degs = [M[centre_indices[0], i] for i in centre_indices]
                level = LCM(core_degs)
                if level.is_squarefree():
                    r = len(level.prime_divisors())
                    cert = {'CM': ZZ(0), 'core_poly':f, 'rho':rho, 'r':r, 'N':level, 'core_degs':core_degs}
                    return True, cert
            print("No central curve found")
        else:
            return True

    # Now we are undecided.  This can happen if either (1) E is not a
    # Q-curve but we did not use enough primes in Step 4 to detect
    # this, or (2) E is a Q-curve but in Step 5 we did not compute the
    # complete isogeny class.  Case (2) is most unlikely since the
    # heuristic bound used in computing isogeny classes means that we
    # have all isogenous curves linked to E by an isogeny of degree
    # supported on primes<1000.

    # We first rerun Step 4 with a larger bound.

    xmaxp = 10*maxp
    if verbose:
        print("Undecided after first round, so we apply more local tests, up to {}".format(xmaxp))

    res4, p = Step4Test(E, B=xmaxp, oldB=maxp, verbose=verbose)
    if not res4:
        if verbose:
            print("No: local test at p={} failed".format(p))
        if certificate:
            return False, None
        else:
            return False

    # Now we rerun Step 5 using a rigorous computaion of the complete
    # isogeny class.  This will probably contain no more curves than
    # before, in which case -- since we already tested that the set of
    # j-invariants does not contain a complete Galois conjugacy class
    # -- we can deduce that E is not a Q-curve.

    if verbose:
        print("...all local tests pass for p<={}".format(xmaxp))
        print("We now compute the complete isogeny class...")

    Cfull = E.isogeny_class(minimal_models=False)
    jCfull = [E2.j_invariant() for E2 in Cfull]

    if len(jC) == len(jCfull):
        if verbose:
            print("...and find that we already had the complete class:so No")
        if certificate:
            return False, None
        else:
            return False
    if verbose:
        print("...and find that the class contains {} curves, not just the {} we computed originally".format(len(jCfull), len(jC)))
    res, cert = conjugacy_test(jCfull, verbose=verbose)
    if res:
        if verbose:
            print("Yes: the isogeny class contains a complete conjugacy class of j-invariants")
        if certificate:
            return True, cert
        else:
            return True
    if verbose:
        print("No: the isogeny class does *not* contain a complete conjugacy class of j-invariants")
    if certificate:
        return False, None
    else:
        return False

def Step4Test(E, B, oldB=0, verbose=False):
    K = E.base_field()
    NN = E.conductor().norm()
    for p in primes(B):
        if p<=oldB or p.divides(NN):
            continue
        Plist = K.primes_above(p)
        if len(Plist)<2:
            continue

        EmodP = [E.reduction(P) for P in Plist]

        # (a) Check all are ordinary or all supersingular:
        ordinary = [Ei.is_ordinary() for Ei in EmodP]
        consistent = all(ordinary) or not any(ordinary)
        if not consistent:
            if verbose:
                print("No: inconsistency at the {} primes dividing {} ".format(len(Plist),p))
                print("  - ordinary: {}".format(ordinary))
            return False, p

        # (b) Skip if all are supersingular:
        if not ordinary[0]:
            continue

        # else compare a_P^2-4*N(P) which should have the same squarefree part:

        discs = [(Ei.trace_of_frobenius()**2 - 4*P.norm()).squarefree_part() for P,Ei in zip(Plist, EmodP)]
        if any([d != discs[0] for d in discs[1:]]):
            if verbose:
                print("No: inconsistency at the {} ordinary primes dividing {} ".format(len(Plist),p))
                print("  - Frobenius discrimants mod squares: {}".format(discs))
            return False, p
    # Now we have failed to prove that E is not a Q-curve
    return True, 0

def conjugacy_test(jC, verbose=True):
    r"""Return True if a list of j-invariants contains a complete
    conjugacy class of 2-power degree.

    - `jC` (list): a list of algebraic numbers in the same field

    Output is (``False``, ``None``) or (``True``, poly) where poly is
    an irreducible polynomial of 2-power degree all of whose roots are
    in the list.
    """
    jQ = next((j for j in jC if j in QQ), None)
    if jQ:
        if verbose:
            print("Yes: an isogenous curve has rational j-invariant {}".format(jQ))
        x = polygen(QQ)
        return True, [x-jQ]

    # If the degree d is odd then we know that none of the
    # j-invariants in the class have 2-power degree, so we can exit.

    K = jC[0].parent()
    if K.degree()%2:
        if verbose:
            print("Odd-degree case: no rational j-invariant in the class {}".format(jC))
        return False, None

    # If K has no quadratic subfields we can similarly conclude right
    # away.  This is one way of determining this.

    if K(1).descend_mod_power(QQ,2) == [1]:
        if verbose:
            print("No-quadratic-subfield case: no rational j-invariant in the class {}".format(jC))
        return False, None

    # compute the minimum polynomials of the j-invariants in the class
    pols = [j.minpoly() for j in jC]

    # pick out those of 2-power degree
    pols = [f for f in pols if f.degree().prime_to_m_part(2)==1]

    # see if there's a poly of degree d appearing d times.  NB There
    # may be more than one of these, possibly including some conjugacy
    # classes defined over the core field but not central, so we
    # return all those with the minimal degree.

    mindeg = min([f.degree() for f in pols])
    minpols = [f for f in pols if f.degree()==mindeg]
    centrepols = list(Set([f for f in pols if f.degree()==minpols.count(f)]))
    if centrepols:
        if verbose:
            print("Yes: the isogeny class contains all j-invariants with min poly {}".format(centrepols))
        return True, centrepols
    if verbose:
        print("No complete conjugacy class of 2-power size found in {}".format(jC))
    return False, None
