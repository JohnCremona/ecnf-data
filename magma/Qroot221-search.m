print "Field 2.2.221.1";
Qx<x> := PolynomialRing(RationalField());
K<w> := NumberField(x^2 - x - 55);
OK := Integers(K);
Plist := [];
Append(~Plist,(2)*OK);
Append(~Plist,(5)*OK+(w)*OK);
Append(~Plist,(5)*OK+(w + 4)*OK);
Append(~Plist,(7)*OK+(w + 2)*OK);
Append(~Plist,(7)*OK+(w + 4)*OK);
Append(~Plist,(3)*OK);
Append(~Plist,(11)*OK+(w)*OK);
Append(~Plist,(11)*OK+(w + 10)*OK);
Append(~Plist,(w + 6)*OK);
Append(~Plist,(w + 8)*OK);
Append(~Plist,(31)*OK+(w + 14)*OK);
Append(~Plist,(31)*OK+(w + 16)*OK);
Append(~Plist,(37)*OK+(w + 15)*OK);
Append(~Plist,(37)*OK+(w + 21)*OK);
Append(~Plist,(41)*OK+(w + 18)*OK);
Append(~Plist,(41)*OK+(w + 22)*OK);
Append(~Plist,(w + 3)*OK);
Append(~Plist,(w - 4)*OK);
Append(~Plist,(w + 1)*OK);
Append(~Plist,(w - 2)*OK);
Append(~Plist,(71)*OK+(w + 23)*OK);
Append(~Plist,(71)*OK+(w + 47)*OK);
Append(~Plist,(73)*OK+(w + 20)*OK);
Append(~Plist,(73)*OK+(w + 52)*OK);
Append(~Plist,(97)*OK+(w + 33)*OK);
Append(~Plist,(97)*OK+(w + 63)*OK);
Append(~Plist,(w + 12)*OK);
Append(~Plist,(w - 13)*OK);
Append(~Plist,(-2*w + 19)*OK);
Append(~Plist,(2*w + 17)*OK);
Append(~Plist,(109)*OK+(w + 24)*OK);
Append(~Plist,(109)*OK+(w + 84)*OK);
Append(~Plist,(w + 13)*OK);
Append(~Plist,(w - 14)*OK);
Append(~Plist,(-2*w + 9)*OK);
Append(~Plist,(2*w + 7)*OK);
Append(~Plist,(163)*OK+(w + 29)*OK);
Append(~Plist,(163)*OK+(w + 133)*OK);
Append(~Plist,(167)*OK+(w + 43)*OK);
Append(~Plist,(167)*OK+(w + 123)*OK);
Append(~Plist,(-2*w + 21)*OK);
Append(~Plist,(2*w + 19)*OK);
Append(~Plist,(-3*w - 16)*OK);
Append(~Plist,(3*w - 19)*OK);
Append(~Plist,(193)*OK+(w + 37)*OK);
Append(~Plist,(193)*OK+(w + 155)*OK);
Append(~Plist,(197)*OK+(w + 78)*OK);
Append(~Plist,(197)*OK+(w + 118)*OK);
Append(~Plist,(227)*OK+(w + 34)*OK);
Append(~Plist,(227)*OK+(w + 192)*OK);
Append(~Plist,(241)*OK+(w + 35)*OK);
Append(~Plist,(241)*OK+(w + 205)*OK);
Append(~Plist,(w + 17)*OK);
Append(~Plist,(w - 18)*OK);
Append(~Plist,(-3*w + 17)*OK);
Append(~Plist,(3*w + 14)*OK);
Append(~Plist,(-2*w + 23)*OK);
Append(~Plist,(2*w + 21)*OK);
Append(~Plist,(317)*OK+(w + 40)*OK);
Append(~Plist,(317)*OK+(w + 276)*OK);
Append(~Plist,(19)*OK);
Append(~Plist,(-3*w + 31)*OK);
Append(~Plist,(3*w + 28)*OK);
Append(~Plist,(379)*OK+(w + 166)*OK);
Append(~Plist,(379)*OK+(w + 212)*OK);
Append(~Plist,(-5*w + 34)*OK);
Append(~Plist,(5*w + 29)*OK);
Append(~Plist,(397)*OK+(w + 66)*OK);
Append(~Plist,(397)*OK+(w + 330)*OK);
Append(~Plist,(401)*OK+(w + 53)*OK);
Append(~Plist,(401)*OK+(w + 347)*OK);
Append(~Plist,(431)*OK+(w + 126)*OK);
Append(~Plist,(431)*OK+(w + 304)*OK);
Append(~Plist,(-3*w - 29)*OK);
Append(~Plist,(3*w - 32)*OK);
Append(~Plist,(-4*w + 23)*OK);
Append(~Plist,(4*w + 19)*OK);
Append(~Plist,(449)*OK+(w + 195)*OK);
Append(~Plist,(449)*OK+(w + 253)*OK);
Append(~Plist,(-3*w - 4)*OK);
Append(~Plist,(3*w - 7)*OK);
Append(~Plist,(479)*OK+(w + 49)*OK);
Append(~Plist,(479)*OK+(w + 429)*OK);
Append(~Plist,(487)*OK+(w + 141)*OK);
Append(~Plist,(487)*OK+(w + 345)*OK);
Append(~Plist,(-3*w - 1)*OK);
Append(~Plist,(3*w - 4)*OK);
Append(~Plist,(499)*OK+(w + 50)*OK);
Append(~Plist,(499)*OK+(w + 448)*OK);
Append(~Plist,(-4*w - 17)*OK);
Append(~Plist,(4*w - 21)*OK);
Append(~Plist,(23)*OK);
Append(~Plist,(541)*OK+(w + 77)*OK);
Append(~Plist,(541)*OK+(w + 463)*OK);
Append(~Plist,(-2*w + 29)*OK);
Append(~Plist,(2*w + 27)*OK);
Append(~Plist,(-5*w - 26)*OK);
Append(~Plist,(5*w - 31)*OK);
Append(~Plist,(5*w + 42)*OK);
Append(~Plist,(5*w - 47)*OK);
effort := 400;
ECSearch := procedure(class_label, N, aplist);
print "Isogeny class ", class_label;
goodP := [P: P in Plist | Valuation(N,P) eq 0];
goodP := [goodP[i]: i in [1..#(aplist)]];
curves := EllipticCurveSearch(N,effort : Primes:=goodP, Traces:=aplist);
curves := [E: E in curves | &and[TraceOfFrobenius(E,goodP[i]) eq aplist[i] : i in [1..#(aplist)]]];
if #curves eq 0 then print "No curve found"; end if;
for E in curves do;
 a1,a2,a3,a4,a6:=Explode(aInvariants(E));
 printf "Curve [%o,%o,%o,%o,%o]\n",a1,a2,a3,a4,a6;
 end for;
end procedure;
SetColumns(0);

ECSearch("4.1-a",(2)*OK,[-3, 3, 3, -3, -5, 0, 0, -4, 3, 0, 0, -3, 3, 0, 0, 1, 1, 6, 6, -15, 15, -6, 6, 12, -12, 12, 12, 14, 14]);
ECSearch("4.1-b",(2)*OK,[3, -3, -3, 3, -5, 0, 0, -4, 3, 0, 0, 3, -3, 0, 0, 1, 1, 6, 6, 15, -15, 6, -6, -12, 12, 12, 12, 14, 14]);
ECSearch("9.1-a",(3)*OK,[-1, -4, 0, 4, 0, -6, 2, 2, -2, 0, 4, 10, 2, 0, 4, 4, -4, -6, -14, 2, 10, 14, -10, -2, 6, 14, -2, -8, -8]);
ECSearch("9.1-b",(3)*OK,[0, -3, 3, -2, 2, 5, -5, 1, 8, 10, -10, 2, -2, 5, -5, 1, 1, 6, 6, 0, 0, -6, 6, -8, 8, 12, 12, 9, 9]);
ECSearch("9.1-c",(3)*OK,[-1, 0, 4, 0, -4, -2, 6, 2, -2, -4, 0, -2, -10, -4, 0, -4, 4, -14, -6, -10, -2, 10, -14, -6, 2, -2, 14, -8, -8]);
ECSearch("9.1-d",(3)*OK,[-3, 0, 0, 4, -4, -4, 4, -2, 2, 4, -4, 8, -8, 8, -8, 4, 4, -6, -6, 12, -12, 0, 0, 16, -16, -6, -6, 0, 0]);
ECSearch("9.1-e",(3)*OK,[-3, 0, 0, -4, 4, 4, -4, -2, 2, -4, 4, -8, 8, -8, 8, 4, 4, -6, -6, -12, 12, 0, 0, -16, 16, -6, -6, 0, 0]);
ECSearch("9.1-f",(3)*OK,[-1, 0, -4, 0, 4, 2, -6, 2, -2, 4, 0, 2, 10, 4, 0, -4, 4, -14, -6, 10, 2, -10, 14, 6, -2, -2, 14, -8, -8]);
ECSearch("9.1-g",(3)*OK,[0, 3, -3, 2, -2, -5, 5, 1, 8, -10, 10, -2, 2, -5, 5, 1, 1, 6, 6, 0, 0, 6, -6, 8, -8, 12, 12, 9, 9]);
ECSearch("9.1-h",(3)*OK,[-1, 4, 0, -4, 0, 6, -2, 2, -2, 0, -4, -10, -2, 0, -4, 4, -4, -6, -14, -2, -10, -14, 10, 2, -6, 14, -2, -8, -8]);
ECSearch("13.1-a",(w + 6)*OK,[-3, 2, -2, 0, 0, 2, 0, 0, 2, 4, -4, 10, -10, 2, -2, 4, 4, -2, -2, 0, 0, 10, -10, 14, -14, 6, 6, 0, 0]);
ECSearch("13.1-b",(w + 6)*OK,[-3, -2, 2, 0, 0, 2, 0, 0, 2, -4, 4, -10, 10, -2, 2, 4, 4, -2, -2, 0, 0, -10, 10, -14, 14, 6, 6, 0, 0]);
ECSearch("17.1-a",(w + 8)*OK,[1, 2, -2, 4, 0, 2, 0, -4, -2, 0, -4, -2, 2, -6, -2, -4, 4, 6, -10, 12, 0, 2, -10, -6, 14, 6, -10, -16, 0]);
ECSearch("17.1-b",(w + 8)*OK,[1, 2, -2, 0, -4, 2, 4, 0, -2, 4, 0, -2, 2, 2, 6, 4, -4, -10, 6, 0, -12, 10, -2, -14, 6, -10, 6, 0, -16]);
ECSearch("17.1-c",(w + 8)*OK,[-3, 2, 2, -4, -4, -6, 0, 0, -2, -4, -4, 2, 2, 6, 6, 4, 4, 6, 6, 4, 4, 6, 6, -2, -2, -10, -10, 8, 8]);
ECSearch("17.1-d",(w + 8)*OK,[-3, -2, -2, 4, 4, -6, 0, 0, -2, 4, 4, -2, -2, -6, -6, 4, 4, 6, 6, -4, -4, -6, -6, 2, 2, -10, -10, 8, 8]);
ECSearch("17.1-e",(w + 8)*OK,[1, -2, 2, 0, 4, 2, -4, 0, -2, -4, 0, 2, -2, -2, -6, 4, -4, -10, 6, 0, 12, -10, 2, 14, -6, -10, 6, 0, -16]);
ECSearch("17.1-f",(w + 8)*OK,[1, -2, 2, -4, 0, 2, 0, 4, -2, 0, 4, 2, -2, 6, 2, -4, 4, 6, -10, -12, 0, -2, 10, 6, -14, 6, -10, -16, 0]);