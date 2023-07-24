r"""A function to check whether an elliptic curve defined over a
totally real field is modular, based on the Magma code by Samir Siksek
in modularitycheck.m.  That file states:

This is based on theorems in the following papers

[Chen] Imin Chen, The Jacobian of Modular Curves Associated to Cartan
Subgroups, Oxford DPhil thesis, 1996.

[FLS] Nuno Freitas, Bao Le Hung, and S. Siksek, Elliptic curves over
Real Quadratic Fields are Modular, Invent. math. (2015) 201,
pp. 159--206.

[Thorne] Jack Thorne, Automorphy of some residually dihedral Galois
representations.  Mathematische Annalen 364 (2016), No. 1--2,
pp. 589--648

Given an elliptic curve $E$ over a totally real field $K$ this returns
either True or False.  True means that E is proven to be modular.
False means that E might not yet be proven to be modular.

The function does nothing for curves over fields K which are not
totally real, returning False.  Otherwise, it returns True for fields
of degree 1 and 2 (with no work needed); else it checks that the
j-invariant does not satisfy a number of properties.  Any curve which
fails this has the (interesting) propery that its mod-p Galois
representations for p=3, 5 and 7 all have certain images (which has
been proved to be impossible for degrees less than 3, as the
associated modular curves have no points of degree < 3 other than
cusps and CM points.

In more detail, E is modular unless for all three of p=3, p=5 and p=7,
j(E) is the image of a K-rational point on at least one of the modular
curves X_0(p), X_s(p), X_ns(p).

"""
from sage.all import (QQ, polygen)
from sage.schemes.elliptic_curves.isogeny_small_degree import Fricke_polynomial

def modular_polynomial(group, p, j):
    x = polygen(j.parent())
    assert (group, p) in [("borel", 3), ("borel", 5), ("borel", 7),
                          ("split", 3), ("split", 5), ("split", 7),
                          ("nonsplit", 5), ("nonsplit", 7)]
    if group == "borel":
        return Fricke_polynomial(p)(x)-j*x
    if group == "split":
        if p == 3:
            return (x-9)**3*(x+3)**3-j*x**3
        if p == 5:
            return ((x**2-5)*(x**2+5*x+10)*(x+5))**3-j*(x**2+5*x+5)**5
        if p == 7:
            return ((x**2-5*x+8)*(x**2-5*x+1)*(x**4-5*x**3+8*x**2-7*x+7)*(x+1))**3*x-j*(x**3-4*x**2+3*x+1)**7
    if group == "nonsplit":
        if p == 5:
            return 5**4*(2*x+1)*(x+1)**3*(6*x**2+21*x+19)**3-j*(x**2+x-1)**5
        if p == 7:
            return ((4*x**2+5*x+2)*(x**2+3*x+4)*(x**2+10*x+4)*(3*x+1))**3-j*(x**3+x**2-2*x-1)**7
    raise RuntimeError # cannot happen -- the assertion would have failed.

def not_from(group, p, j):
    return len(modular_polynomial(group, p, j).roots()) == 0

def isModular(E):
    K = E.base_field()
    if K is QQ:
        return True
    if not K.is_totally_real():
        return False # unknown
    if K.degree() == 2:
        return True
    if E.has_cm():
        return True

    j = E.j_invariant()

    if not_from("borel", 3, j) and not_from("split", 3, j):
        return True

    if not_from("borel", 5, j) and (not K(5).is_square() or (not_from("split", 5, j) and not_from("nonsplit", 5, j))):
        return True

    if not_from("borel", 7, j) and not_from("split", 7, j) and not_from("nonsplit", 7, j):
        return True

    return False # We've run out of tricks!
