B=10000

print "Qsqrt-1"

E1f1 = list(curves_K1(B,1))
E1f2 = list(curves_K1(B,2))
E1 = E1f1 + E1f2
assert all_non_iso(E1)
print "   found %s+%s=%s curves" % (len(E1f1),len(E1f2),len(E1))

print "Qsqrt-2"

E2 = E2f1 = list(curves_K2(B))
assert all_non_iso(E2)
print "   found %s curves" % len(E2)

print "Qsqrt-3"

E3f1 = list(curves_K3(B,1))
E3f2 = list(curves_K3(B,2))
E3f3 = list(curves_K3(B,3))
E3 = E3f1+E3f2+E3f3
assert all_non_iso(E3)
print "   found %s+%s+%s=%s curves" % (len(E3f1),len(E3f2),len(E3f3),len(E3))

print "Qsqrt-7"

E7f1 = list(curves_K7(B,1))
assert all_non_iso(E7f1)
E7f2 = list(curves_K7(B,2))
assert all_non_iso(E7f2)
E7 = E7f1+E7f2
print "   found %s+%s=%s curves" % (len(E7f1),len(E7f2),len(E7))

print "Qsqrt-11"

E11 = E11f1 = list(curves_K11(B))
assert all_non_iso(E11)
print "   found %s curves" % len(E11)
