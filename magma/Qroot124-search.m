print "Field 2.2.124.1";
Qx<x> := PolynomialRing(RationalField());
K<w> := NumberField(x^2 - 31);
OK := Integers(K);
Plist := [];
Append(~Plist,(7*w + 39)*OK);
Append(~Plist,(2*w + 11)*OK);
Append(~Plist,(2*w - 11)*OK);
Append(~Plist,(-w - 6)*OK);
Append(~Plist,(w - 6)*OK);
Append(~Plist,(9*w - 50)*OK);
Append(~Plist,(-9*w - 50)*OK);
Append(~Plist,(3*w - 16)*OK);
Append(~Plist,(-3*w - 16)*OK);
Append(~Plist,(w)*OK);
Append(~Plist,(8*w - 45)*OK);
Append(~Plist,(8*w + 45)*OK);
Append(~Plist,(-2*w + 9)*OK);
Append(~Plist,(2*w + 9)*OK);
Append(~Plist,(7)*OK);
Append(~Plist,(-20*w + 111)*OK);
Append(~Plist,(-20*w - 111)*OK);
Append(~Plist,(-3*w + 14)*OK);
Append(~Plist,(3*w + 14)*OK);
Append(~Plist,(-33*w + 184)*OK);
Append(~Plist,(33*w + 184)*OK);
Append(~Plist,(-2*w + 15)*OK);
Append(~Plist,(2*w + 15)*OK);
Append(~Plist,(-6*w - 35)*OK);
Append(~Plist,(-6*w + 35)*OK);
Append(~Plist,(w + 12)*OK);
Append(~Plist,(w - 12)*OK);
Append(~Plist,(-41*w + 228)*OK);
Append(~Plist,(-41*w - 228)*OK);
Append(~Plist,(-62*w + 345)*OK);
Append(~Plist,(55*w - 306)*OK);
Append(~Plist,(10*w - 57)*OK);
Append(~Plist,(10*w + 57)*OK);
Append(~Plist,(11*w - 60)*OK);
Append(~Plist,(-11*w - 60)*OK);
Append(~Plist,(18*w - 101)*OK);
Append(~Plist,(18*w + 101)*OK);
Append(~Plist,(24*w - 133)*OK);
Append(~Plist,(24*w + 133)*OK);
Append(~Plist,(13)*OK);
Append(~Plist,(29*w - 162)*OK);
Append(~Plist,(29*w + 162)*OK);
Append(~Plist,(-3*w - 10)*OK);
Append(~Plist,(3*w - 10)*OK);
Append(~Plist,(-9*w + 52)*OK);
Append(~Plist,(9*w + 52)*OK);
Append(~Plist,(5*w - 24)*OK);
Append(~Plist,(5*w + 24)*OK);
Append(~Plist,(7*w + 36)*OK);
Append(~Plist,(-7*w + 36)*OK);
Append(~Plist,(-4*w - 27)*OK);
Append(~Plist,(4*w - 27)*OK);
Append(~Plist,(12*w + 65)*OK);
Append(~Plist,(12*w - 65)*OK);
Append(~Plist,(15*w - 82)*OK);
Append(~Plist,(-15*w - 82)*OK);
Append(~Plist,(-17*w - 96)*OK);
Append(~Plist,(-17*w + 96)*OK);
Append(~Plist,(-3*w - 4)*OK);
Append(~Plist,(3*w - 4)*OK);
Append(~Plist,(-4*w + 15)*OK);
Append(~Plist,(4*w + 15)*OK);
Append(~Plist,(-43*w - 240)*OK);
Append(~Plist,(-43*w + 240)*OK);
Append(~Plist,(17)*OK);
Append(~Plist,(w + 18)*OK);
Append(~Plist,(w - 18)*OK);
Append(~Plist,(-2*w + 21)*OK);
Append(~Plist,(2*w + 21)*OK);
Append(~Plist,(25*w + 138)*OK);
Append(~Plist,(-25*w + 138)*OK);
Append(~Plist,(66*w + 367)*OK);
Append(~Plist,(-66*w + 367)*OK);
Append(~Plist,(-75*w + 418)*OK);
Append(~Plist,(75*w + 418)*OK);
Append(~Plist,(19)*OK);
Append(~Plist,(16*w + 87)*OK);
Append(~Plist,(16*w - 87)*OK);
Append(~Plist,(-39*w + 218)*OK);
Append(~Plist,(39*w + 218)*OK);
Append(~Plist,(87*w + 484)*OK);
Append(~Plist,(87*w - 484)*OK);
Append(~Plist,(-3*w - 26)*OK);
Append(~Plist,(3*w - 26)*OK);
Append(~Plist,(-15*w - 86)*OK);
Append(~Plist,(15*w - 86)*OK);
Append(~Plist,(8*w - 39)*OK);
Append(~Plist,(8*w + 39)*OK);
Append(~Plist,(-4*w + 3)*OK);
Append(~Plist,(4*w + 3)*OK);
Append(~Plist,(-6*w - 25)*OK);
Append(~Plist,(6*w - 25)*OK);
Append(~Plist,(10*w + 51)*OK);
Append(~Plist,(-10*w + 51)*OK);
Append(~Plist,(-5*w + 36)*OK);
Append(~Plist,(5*w + 36)*OK);
Append(~Plist,(23*w - 126)*OK);
Append(~Plist,(-23*w - 126)*OK);
Append(~Plist,(-42*w - 235)*OK);
Append(~Plist,(-42*w + 235)*OK);
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

ECSearch("3.1-a",(2*w + 11)*OK,[1, 1, -2, 3, 1, 6, 3, 8, -1, 2, 7, 6, -4, -11, 1, 11, 15, -10, 4, -6, -3, 12, -1, 9, 6, -4, -4, -14, -2]);
ECSearch("3.1-b",(2*w + 11)*OK,[1, -1, -2, 3, -1, -6, -3, -8, 1, 2, 7, -6, 4, -11, -1, -11, -15, 10, 4, -6, -3, 12, -1, 9, 6, -4, 4, 14, 2]);
ECSearch("3.2-a",(2*w - 11)*OK,[1, 1, 3, -2, 6, 1, 8, 3, -1, 7, 2, -4, 6, -11, 11, 1, -10, 15, -6, 4, 12, -3, 9, -1, -4, 6, -14, -4, -12]);
ECSearch("3.2-b",(2*w - 11)*OK,[1, -1, 3, -2, -6, -1, -8, -3, 1, 7, 2, 4, -6, -11, -11, -1, 10, -15, -6, 4, 12, -3, 9, -1, -4, 6, 14, 4, 12]);
ECSearch("4.1-a",(2)*OK,[-2, 2, -1, -1, 4, -4, 6, -6, 0, -5, -5, -10, 10, -13, -10, 10, -4, 4, 7, 7, 7, 7, -9, -9, -13, -13, 0, 0, 8]);
ECSearch("4.1-b",(2)*OK,[2, -2, -1, -1, -4, 4, -6, 6, 0, -5, -5, 10, -10, -13, 10, -10, 4, -4, 7, 7, 7, 7, -9, -9, -13, -13, 0, 0, -8]);
ECSearch("5.1-a",(-w - 6)*OK,[-1, -2, 2, -2, -2, -2, 8, 4, 4, 6, -10, 2, -10, 2, -16, 0, 6, -6, -18, -10, -18, -10, -10, -2, 6, -2, -4, -16, -2]);
ECSearch("5.1-b",(-w - 6)*OK,[-1, 2, -2, -2, 2, 2, -8, -4, -4, 6, -10, -2, 10, 2, 16, 0, -6, 6, -18, -10, -18, -10, -10, -2, 6, -2, 4, 16, 2]);
ECSearch("5.2-a",(w - 6)*OK,[-1, 2, -2, -2, -2, -2, 4, 8, 4, -10, 6, -10, 2, 2, 0, -16, -6, 6, -10, -18, -10, -18, -2, -10, -2, 6, -16, -4, 6]);
ECSearch("5.2-b",(w - 6)*OK,[-1, -2, 2, -2, 2, 2, -4, -8, -4, -10, 6, 10, -2, 2, 0, 16, 6, -6, -10, -18, -10, -18, -2, -10, -2, 6, 16, 4, -6]);
ECSearch("6.1-a",(-w + 5)*OK,[2, 2, -2, 0, 2, 6, -2, -8, 2, -6, 10, 0, -10, 10, -10, 6, 8, -10, 2, 6, 6, -10, 14, -6, 2, 6, 2, 0]);
ECSearch("6.1-b",(-w + 5)*OK,[1, 2, -1, 5, -2, 3, 0, 3, 10, -5, 10, 8, 5, 1, 11, 7, 14, 4, -2, -11, -8, 3, 9, 14, 8, 0, 6, -14]);
ECSearch("6.1-c",(-w + 5)*OK,[-1, 2, -1, -5, 2, -3, 0, -3, 10, -5, -10, -8, 5, -1, -11, -7, -14, 4, -2, -11, -8, 3, 9, 14, 8, 0, -6, 14]);
ECSearch("6.1-d",(-w + 5)*OK,[-2, 2, -2, 0, -2, -6, 2, 8, 2, -6, -10, 0, -10, -10, 10, -6, -8, -10, 2, 6, 6, -10, 14, -6, 2, -6, -2, 0]);
ECSearch("6.1-e",(-w + 5)*OK,[1, -4, 1, 3, -2, -9, -4, -7, 2, -3, -4, 6, 5, 5, -5, -9, -14, 8, 8, -3, -18, 5, 5, -6, 14, -18, -8, 0]);
ECSearch("6.1-f",(-w + 5)*OK,[-1, -4, 1, -3, 2, 9, 4, 7, 2, -3, 4, -6, 5, -5, 5, 9, 14, 8, 8, -3, -18, 5, 5, -6, 14, 18, 8, 0]);
ECSearch("6.2-a",(w + 5)*OK,[2, -2, 2, 2, 0, -2, 6, -8, -6, 2, 0, 10, -10, -10, 10, 8, 6, 2, -10, 6, 6, 14, -10, 2, -6, 2, 6, -14]);
ECSearch("6.2-b",(w + 5)*OK,[1, -1, 2, -2, 5, 0, 3, 3, -5, 10, 8, 10, 5, 11, 1, 14, 7, -2, 4, -8, -11, 9, 3, 8, 14, 6, 0, 8]);
ECSearch("6.2-c",(w + 5)*OK,[-1, -1, 2, 2, -5, 0, -3, -3, -5, 10, -8, -10, 5, -11, -1, -14, -7, -2, 4, -8, -11, 9, 3, 8, 14, -6, 0, -8]);
ECSearch("6.2-d",(w + 5)*OK,[-2, -2, 2, -2, 0, 2, -6, 8, -6, 2, 0, -10, -10, 10, -10, -8, -6, 2, -10, 6, 6, 14, -10, 2, -6, -2, -6, 14]);
ECSearch("6.2-e",(w + 5)*OK,[1, 1, -4, -2, 3, -4, -9, -7, -3, 2, 6, -4, 5, -5, 5, -14, -9, 8, 8, -18, -3, 5, 5, 14, -6, -8, -18, -10]);
ECSearch("6.2-f",(w + 5)*OK,[-1, 1, -4, 2, -3, 4, 9, 7, -3, 2, -6, 4, 5, 5, -5, 14, 9, 8, 8, -18, -3, 5, 5, 14, -6, 8, 18, 10]);
ECSearch("8.1-a",(14*w + 78)*OK,[-2, 2, 3, 3, 4, -4, -2, 2, 0, 11, 11, 6, -6, -13, 14, -14, 12, -12, -1, -1, 3, 3, 3, 3, 11, 11, 0, 0, -8]);
ECSearch("8.1-b",(14*w + 78)*OK,[2, -2, 3, 3, -4, 4, 2, -2, 0, 11, 11, -6, 6, -13, -14, 14, -12, 12, -1, -1, 3, 3, 3, 3, 11, 11, 0, 0, 8]);
ECSearch("8.1-c",(14*w + 78)*OK,[-2, 2, -2, -2, -6, 6, 8, -8, 0, 6, 6, 6, -6, 2, 4, -4, 2, -2, 14, 14, -2, -2, -2, -2, 6, 6, 20, -20, 2]);
ECSearch("8.1-d",(14*w + 78)*OK,[2, -2, -2, -2, 6, -6, -8, 8, 0, 6, 6, -6, 6, 2, -4, 4, -2, 2, 14, 14, -2, -2, -2, -2, 6, 6, -20, 20, -2]);
ECSearch("9.2-a",(-5*w - 28)*OK,[1, 3, 1, 4, -6, 3, -6, 3, 3, -5, 4, -6, -6, 5, 15, -3, -12, -3, -8, -8, 8, -1, 11, -7, 8, 2, 6, 6, -18]);
ECSearch("9.2-b",(-5*w - 28)*OK,[1, -3, 1, 4, 6, -3, 6, -3, -3, -5, 4, 6, 6, 5, -15, 3, 12, 3, -8, -8, 8, -1, 11, -7, 8, 2, -6, -6, 18]);
ECSearch("9.2-c",(-5*w - 28)*OK,[-2, 0, 4, -2, 0, 0, 0, 0, 0, 10, -8, 0, 0, 14, 0, 0, 0, 0, -8, -8, 2, 20, 20, 20, -16, 14, 0, 0, 0]);
ECSearch("9.3-a",(5*w - 28)*OK,[1, 3, 4, 1, 3, -6, 3, -6, 3, 4, -5, -6, -6, 5, -3, 15, -3, -12, -8, -8, -1, 8, -7, 11, 2, 8, 6, 6, 0]);
ECSearch("9.3-b",(5*w - 28)*OK,[1, -3, 4, 1, -3, 6, -3, 6, -3, 4, -5, 6, 6, 5, 3, -15, 3, 12, -8, -8, -1, 8, -7, 11, 2, 8, -6, -6, 0]);
ECSearch("9.3-c",(5*w - 28)*OK,[-2, 0, -2, 4, 0, 0, 0, 0, 0, -8, 10, 0, 0, 14, 0, 0, 0, 0, -8, -8, 20, 2, 20, 20, 14, -16, 0, 0, 0]);
ECSearch("10.1-a",(-3*w + 17)*OK,[0, 2, -2, -4, 2, -6, 2, -8, -6, -6, 2, -4, 6, 6, 10, 6, 12, -10, 2, 6, -18, -2, -18, -6, 2, 10, 22, -4]);
ECSearch("10.1-b",(-3*w + 17)*OK,[0, -2, -2, 4, -2, 6, -2, 8, -6, -6, -2, 4, 6, -6, -10, -6, -12, -10, 2, 6, -18, -2, -18, -6, 2, -10, -22, 4]);
ECSearch("10.2-a",(-3*w - 17)*OK,[2, 0, -2, 2, -4, 2, -6, -8, -6, -6, -4, 2, 6, 10, 6, 12, 6, 2, -10, -18, 6, -18, -2, 2, -6, 22, 10, 2]);
ECSearch("10.2-b",(-3*w - 17)*OK,[-2, 0, -2, -2, 4, -2, 6, 8, -6, -6, 4, -2, 6, -10, -6, -12, -6, 2, -10, -18, 6, -18, -2, 2, -6, -22, -10, -2]);
ECSearch("12.1-a",(4*w + 22)*OK,[-1, 2, 1, 3, 2, 3, 4, -5, -10, 3, -2, 12, 5, 1, 11, 9, 2, 8, 2, 3, 0, 5, -1, 6, 8, -12, 14, -6]);
ECSearch("12.1-b",(4*w + 22)*OK,[1, 2, 1, -3, -2, -3, -4, 5, -10, 3, 2, -12, 5, -1, -11, -9, -2, 8, 2, 3, 0, 5, -1, 6, 8, 12, -14, 6]);
ECSearch("12.2-a",(4*w - 22)*OK,[-1, 1, 2, 2, 3, 4, 3, -5, 3, -10, 12, -2, 5, 11, 1, 2, 9, 2, 8, 0, 3, -1, 5, 8, 6, 14, -12, 4]);
ECSearch("12.2-b",(4*w - 22)*OK,[1, 1, 2, -2, -3, -4, -3, 5, 3, -10, -12, 2, 5, -11, -1, -2, -9, 2, 8, 0, 3, -1, 5, 8, 6, -14, 12, -4]);
ECSearch("15.2-a",(-w - 4)*OK,[-2, 1, 1, -2, 3, 6, -4, -7, 2, -8, -9, -4, 10, 10, -10, 6, -4, 13, -12, -3, 12, -10, 0, -6, -1, -13, -8, -20]);
ECSearch("15.2-b",(-w - 4)*OK,[-2, -1, 1, 2, -3, -6, 4, 7, 2, -8, 9, 4, 10, -10, 10, -6, 4, 13, -12, -3, 12, -10, 0, -6, -1, 13, 8, 20]);
ECSearch("15.2-c",(-w - 4)*OK,[2, 1, 1, -2, -5, -2, -4, 1, 2, 0, -1, -4, 2, -6, 14, 6, -12, 5, -12, 5, -12, -2, 8, 10, -9, -13, 8, 12]);
ECSearch("15.2-d",(-w - 4)*OK,[2, -1, 1, 2, 5, 2, 4, -1, 2, 0, 1, 4, 2, 6, -14, -6, 12, 5, -12, 5, -12, -2, 8, 10, -9, 13, -8, -12]);
ECSearch("15.3-a",(-w + 4)*OK,[-2, 1, 1, 3, -2, -4, 6, -7, -8, 2, -4, -9, 10, -10, 10, -4, 6, -12, 13, 12, -3, 0, -10, -1, -6, -8, -13, 15]);
ECSearch("15.3-b",(-w + 4)*OK,[-2, -1, 1, -3, 2, 4, -6, 7, -8, 2, 4, 9, 10, 10, -10, 4, -6, -12, 13, 12, -3, 0, -10, -1, -6, 8, 13, -15]);
ECSearch("15.3-c",(-w + 4)*OK,[2, 1, 1, -5, -2, -4, -2, 1, 0, 2, -4, -1, 2, 14, -6, -12, 6, -12, 5, -12, 5, 8, -2, -9, 10, 8, -13, 7]);
ECSearch("15.3-d",(-w + 4)*OK,[2, -1, 1, 5, 2, 4, 2, -1, 0, 2, 4, 1, 2, -14, 6, 12, -6, -12, 5, -12, 5, 8, -2, -9, 10, -8, 13, -7]);
ECSearch("16.1-a",(4)*OK,[2, 2, 3, 3, 0, 0, -6, -6, -4, 3, 3, 2, 2, 11, -10, -10, 12, 12, -17, -17, 3, 3, -5, -5, 3, 3, 4, 4, -16]);
ECSearch("16.1-b",(4)*OK,[-2, -2, 3, 3, 0, 0, 6, 6, 4, 3, 3, -2, -2, 11, 10, 10, -12, -12, -17, -17, 3, 3, -5, -5, 3, 3, -4, -4, 16]);
ECSearch("18.1-a",(21*w + 117)*OK,[-3, 4, -5, 2, 7, 0, -3, 5, 5, -5, 2, -5, 14, 14, -2, -9, 19, -2, 19, 12, -17, 4, 14, -7, -20, -6, -16]);
ECSearch("18.1-b",(21*w + 117)*OK,[-1, 4, -3, -2, -3, 8, 5, 11, -9, 11, -6, 11, 2, -2, 6, -11, -1, 14, -3, 12, 11, -4, -6, -1, -12, 22, 0]);
ECSearch("18.1-c",(21*w + 117)*OK,[-3, 4, 5, -2, -7, 0, 3, 5, 5, 5, -2, -5, -14, -14, 2, 9, 19, -2, 19, 12, -17, 4, 14, -7, 20, 6, 16]);
ECSearch("18.1-d",(21*w + 117)*OK,[-1, 4, 3, 2, 3, -8, -5, 11, -9, -11, 6, 11, -2, 2, -6, 11, -1, 14, -3, 12, 11, -4, -6, -1, 12, -22, 0]);
ECSearch("18.1-e",(21*w + 117)*OK,[-1, -1, -1, 1, -6, 6, 0, 4, 4, -8, 8, 2, 1, -1, 1, -1, -17, -17, -2, -2, -12, -12, 8, 8, 0, 0, 4]);
ECSearch("18.1-f",(21*w + 117)*OK,[-1, -1, 1, -1, 6, -6, 0, 4, 4, 8, -8, 2, -1, 1, -1, 1, -17, -17, -2, -2, -12, -12, 8, 8, 0, 0, -4]);
ECSearch("18.1-g",(21*w + 117)*OK,[4, -1, -2, -3, 8, -3, 5, -9, 11, -6, 11, 11, -2, 2, -11, 6, 14, -1, 12, -3, -4, 11, -1, -6, 22, -12, 20]);
ECSearch("18.1-h",(21*w + 117)*OK,[4, -3, -2, 5, 0, -7, 3, 5, 5, -2, 5, -5, -14, -14, 9, 2, -2, 19, 12, 19, 4, -17, -7, 14, 6, 20, -12]);
ECSearch("18.1-i",(21*w + 117)*OK,[4, -1, 2, 3, -8, 3, -5, -9, 11, 6, -11, 11, 2, -2, 11, -6, 14, -1, 12, -3, -4, 11, -1, -6, -22, 12, -20]);
ECSearch("18.1-j",(21*w + 117)*OK,[4, -3, 2, -5, 0, 7, -3, 5, 5, 2, -5, -5, 14, 14, -9, -2, -2, 19, 12, 19, 4, -17, -7, 14, -6, -20, 12]);
ECSearch("18.2-a",(-w - 7)*OK,[1, 0, -3, -3, 6, 9, 6, 5, -12, -9, 2, 2, 5, -1, 5, -9, 12, 8, -4, -9, 12, 5, -13, 18, -12, 14, 14, -4]);
ECSearch("18.2-b",(-w - 7)*OK,[1, -2, 3, -3, -2, -1, 4, 7, -2, -5, 2, 4, -3, 13, -9, -9, -6, -4, 2, 1, -8, -9, -11, -18, 8, -4, -18, 14]);
ECSearch("18.2-c",(-w - 7)*OK,[-1, 0, -3, 3, -6, -9, -6, -5, -12, -9, -2, -2, 5, 1, -5, 9, -12, 8, -4, -9, 12, 5, -13, 18, -12, -14, -14, 4]);
ECSearch("18.2-d",(-w - 7)*OK,[-1, -2, 3, 3, 2, 1, -4, -7, -2, -5, -2, -4, -3, -13, 9, 9, 6, -4, 2, 1, -8, -9, -11, -18, 8, 4, 18, -14]);
ECSearch("18.3-a",(-w + 7)*OK,[1, -3, 0, 6, -3, 6, 9, 5, -9, -12, 2, 2, 5, 5, -1, 12, -9, -4, 8, 12, -9, -13, 5, -12, 18, 14, 14, 14]);
ECSearch("18.3-b",(-w + 7)*OK,[1, 3, -2, -2, -3, 4, -1, 7, -5, -2, 4, 2, -3, -9, 13, -6, -9, 2, -4, -8, 1, -11, -9, 8, -18, -18, -4, -4]);
ECSearch("18.3-c",(-w + 7)*OK,[-1, -3, 0, -6, 3, -6, -9, -5, -9, -12, -2, -2, 5, -5, 1, -12, 9, -4, 8, 12, -9, -13, 5, -12, 18, -14, -14, -14]);
ECSearch("18.3-d",(-w + 7)*OK,[-1, 3, -2, 2, 3, -4, 1, -7, -5, -2, -4, -2, -3, 9, -13, 6, 9, 2, -4, -8, 1, -11, -9, 8, -18, 18, 4, 4]);
ECSearch("24.1-a",(-2*w + 10)*OK,[1, 0, -3, -1, -2, -1, 4, -3, 2, -7, 0, -6, 5, 13, -13, 15, -6, -16, -4, -3, 18, 9, -3, 2, 2, -6, 12, 20]);
ECSearch("24.1-b",(-2*w + 10)*OK,[-1, 0, -3, 1, 2, 1, -4, 3, 2, -7, 0, 6, 5, -13, 13, -15, 6, -16, -4, -3, 18, 9, -3, 2, 2, 6, -12, -20]);
ECSearch("24.1-c",(-2*w + 10)*OK,[1, 4, 1, 3, 6, 7, -4, 9, -6, -3, -12, 6, 5, 5, -5, 7, 2, -16, -16, 13, -2, -11, -11, -6, -18, -2, -16, -8]);
ECSearch("24.1-d",(-2*w + 10)*OK,[-1, 4, 1, -3, -6, -7, 4, -9, -6, -3, 12, -6, 5, -5, 5, -7, -2, -16, -16, 13, -2, -11, -11, -6, -18, 2, 16, 8]);
ECSearch("24.2-a",(2*w + 10)*OK,[1, -3, 0, -2, -1, 4, -1, -3, -7, 2, -6, 0, 5, -13, 13, -6, 15, -4, -16, 18, -3, -3, 9, 2, 2, 12, -6, -14]);
ECSearch("24.2-b",(2*w + 10)*OK,[-1, -3, 0, 2, 1, -4, 1, 3, -7, 2, 6, 0, 5, 13, -13, 6, -15, -4, -16, 18, -3, -3, 9, 2, 2, -12, 6, 14]);
ECSearch("24.2-c",(2*w + 10)*OK,[1, 1, 4, 6, 3, -4, 7, 9, -3, -6, 6, -12, 5, -5, 5, 2, 7, -16, -16, -2, 13, -11, -11, -18, -6, -16, -2, -10]);
ECSearch("24.2-d",(2*w + 10)*OK,[-1, 1, 4, -6, -3, 4, -7, -9, -3, -6, -6, 12, 5, 5, -5, -2, -7, -16, -16, -2, 13, -11, -11, -18, -6, 16, 2, 10]);
ECSearch("25.2-a",(-12*w + 67)*OK,[1, 1, 1, -2, -4, 1, 8, -2, 4, 2, -3, 11, 1, -11, 6, 16, 0, 15, 14, -6, 2, -18, -16, 4, 1, -9, -4, 16, -12]);
ECSearch("25.2-b",(-12*w + 67)*OK,[1, -1, -1, -2, 4, -1, -8, 2, -4, 2, -3, -11, -1, -11, -6, -16, 0, -15, 14, -6, 2, -18, -16, 4, 1, -9, 4, -16, 12]);
ECSearch("25.2-c",(-12*w + 67)*OK,[2, 0, 0, 2, 0, 0, 0, 0, 0, -8, -8, 0, 0, 0, 0, 0, 0, 0, -8, -18, 20, -20, 6, 6, -14, 16, 0, 0, 0]);
ECSearch("25.3-a",(12*w + 67)*OK,[1, 1, 1, -2, 1, -4, -2, 8, 4, -3, 2, 1, 11, -11, 16, 6, 15, 0, -6, 14, -18, 2, 4, -16, -9, 1, 16, -4, 13]);
ECSearch("25.3-b",(12*w + 67)*OK,[1, -1, -1, -2, -1, 4, 2, -8, -4, -3, 2, -1, -11, -11, -16, -6, -15, 0, -6, 14, -18, 2, 4, -16, -9, 1, -16, 4, -13]);
ECSearch("25.3-c",(12*w + 67)*OK,[2, 0, 0, 2, 0, 0, 0, 0, 0, -8, -8, 0, 0, 0, 0, 0, 0, 0, -18, -8, -20, 20, 6, 6, 16, -14, 0, 0, 0]);
ECSearch("27.1-a",(6*w + 33)*OK,[1, -3, 1, 3, 1, -4, -6, -4, 4, 2, 8, -6, -2, 11, 1, 11, 9, -15, 13, 6, -6, 0, -10, 2, -6, -8, 20, 0]);
ECSearch("27.1-b",(6*w + 33)*OK,[1, -3, 1, -3, -1, 4, 6, 4, 4, 2, -8, 6, -2, -11, -1, -11, -9, -15, 13, 6, -6, 0, -10, 2, -6, 8, -20, 0]);
ECSearch("27.2-a",(6*w - 33)*OK,[1, 1, -3, 1, 3, -6, -4, -4, 2, 4, -6, 8, -2, 1, 11, 9, 11, 13, -15, -6, 6, -10, 0, -6, 2, 20, -8, 10]);
ECSearch("27.2-b",(6*w - 33)*OK,[1, 1, -3, -1, -3, 6, 4, 4, 2, 4, 6, -8, -2, -1, -11, -9, -11, 13, -15, -6, 6, -10, 0, -6, 2, -20, 8, -10]);
ECSearch("30.1-a",(w + 1)*OK,[1, -3, 3, 0, 3, 6, -7, -6, -9, -10, 8, -1, 11, 5, 9, 0, 8, -10, -9, -6, -7, -7, -18, 12, 8, -22, 20]);
ECSearch("30.1-b",(w + 1)*OK,[-1, -3, -3, 0, -3, -6, 7, -6, -9, 10, -8, -1, -11, -5, -9, 0, 8, -10, -9, -6, -7, -7, -18, 12, -8, 22, -20]);
ECSearch("30.4-a",(-w + 1)*OK,[1, -3, 0, 3, 6, 3, -7, -9, -6, 8, -10, -1, 5, 11, 0, 9, -10, 8, -6, -9, -7, -7, 12, -18, -22, 8, 14]);
ECSearch("30.4-b",(-w + 1)*OK,[-1, -3, 0, -3, -6, -3, 7, -9, -6, -8, 10, -1, -5, -11, 0, -9, -10, 8, -6, -9, -7, -7, 12, -18, 22, -8, -14]);
ECSearch("33.2-a",(-4*w - 23)*OK,[-2, 2, -3, -2, 5, 4, 3, 1, -5, -10, -8, 3, 4, 16, 2, -17, 12, -6, 4, -12, -12, -18, -16, -1, 12, 8, -14, -6]);
ECSearch("33.2-b",(-4*w - 23)*OK,[-2, -2, -3, -2, -5, -4, -3, -1, -5, -10, 8, -3, 4, -16, -2, 17, -12, -6, 4, -12, -12, -18, -16, -1, 12, -8, 14, 6]);
ECSearch("33.3-a",(-4*w + 23)*OK,[-2, 2, -2, -3, 5, 3, 4, 1, -10, -5, 3, -8, 4, 2, 16, 12, -17, 4, -6, -12, -12, -16, -18, 12, -1, -14, 8, -7]);
ECSearch("33.3-b",(-4*w + 23)*OK,[-2, -2, -2, -3, -5, -3, -4, -1, -10, -5, -3, 8, 4, -2, -16, -12, 17, 4, -6, -12, -12, -16, -18, 12, -1, 14, -8, 7]);
ECSearch("36.1-a",(6)*OK,[-3, 1, 2, -2, 2, 6, 0, 5, -7, 2, -2, -5, 4, 4, -4, -12, -17, 7, 13, -15, -5, 19, -19, -7, 2, 6, -14]);
ECSearch("36.1-b",(6)*OK,[2, 2, 4, -4, 0, 0, 0, 10, 10, 8, -8, 2, -16, 16, -4, 4, -2, -2, 10, 10, 18, 18, 2, 2, -12, 12, -16]);
ECSearch("36.1-c",(6)*OK,[1, -3, 2, -2, -6, -2, 0, -7, 5, 2, -2, -5, -4, -4, 12, 4, 7, -17, -15, 13, 19, -5, -7, -19, -6, -2, -6]);
ECSearch("36.1-d",(6)*OK,[-3, 1, -2, 2, -2, -6, 0, 5, -7, -2, 2, -5, -4, -4, 4, 12, -17, 7, 13, -15, -5, 19, -19, -7, -2, -6, 14]);
ECSearch("36.1-e",(6)*OK,[1, -3, -2, 2, 6, 2, 0, -7, 5, -2, 2, -5, 4, 4, -12, -4, 7, -17, -15, 13, 19, -5, -7, -19, 6, 2, 6]);
ECSearch("36.1-f",(6)*OK,[2, 2, -4, 4, 0, 0, 0, 10, 10, -8, 8, 2, 16, -16, 4, -4, -2, -2, 10, 10, 18, 18, 2, 2, 12, -12, 16]);
ECSearch("36.2-a",(-10*w - 56)*OK,[2, -3, 3, 0, 0, -6, 6, 4, 3, -3, -2, -2, 11, 10, 10, -12, 12, 17, 17, 3, -3, 5, 5, -3, 3, 4, 4, 16]);
ECSearch("36.2-b",(-10*w - 56)*OK,[1, -3, -2, -2, -3, -8, -7, 1, -5, 10, -8, 2, -3, -15, -5, -6, -9, 14, -16, -8, -17, -5, -15, -16, 6, 18, 8, 8]);
ECSearch("36.2-c",(-10*w - 56)*OK,[-1, -3, -2, 2, 3, 8, 7, -1, -5, 10, 8, -2, -3, 15, 5, 6, 9, 14, -16, -8, -17, -5, -15, -16, 6, -18, -8, -8]);
ECSearch("36.2-d",(-10*w - 56)*OK,[-2, -3, 3, 0, 0, 6, -6, -4, 3, -3, 2, 2, 11, -10, -10, 12, -12, 17, 17, 3, -3, 5, 5, -3, 3, -4, -4, -16]);
ECSearch("36.2-e",(-10*w - 56)*OK,[1, -1, -4, -2, 5, 6, 3, -9, 7, -8, -10, 10, 5, -1, -11, -4, 7, 16, 4, 4, -11, 9, -9, -4, 2, 18, 6, -22]);
ECSearch("36.2-f",(-10*w - 56)*OK,[-1, -1, -4, 2, -5, -6, -3, 9, 7, -8, 10, -10, 5, 1, 11, 4, -7, 16, 4, 4, -11, 9, -9, -4, 2, -18, -6, 22]);
ECSearch("36.3-a",(10*w - 56)*OK,[2, 3, -3, 0, 0, 6, -6, 4, -3, 3, -2, -2, 11, 10, 10, 12, -12, 17, 17, -3, 3, 5, 5, 3, -3, 4, 4, 16]);
ECSearch("36.3-b",(10*w - 56)*OK,[1, -2, -3, -3, -2, -7, -8, 1, 10, -5, 2, -8, -3, -5, -15, -9, -6, -16, 14, -17, -8, -15, -5, 6, -16, 8, 18, -22]);
ECSearch("36.3-c",(10*w - 56)*OK,[-1, -2, -3, 3, 2, 7, 8, -1, 10, -5, -2, 8, -3, 5, 15, 9, 6, -16, 14, -17, -8, -15, -5, 6, -16, -8, -18, 22]);
ECSearch("36.3-d",(10*w - 56)*OK,[-2, 3, -3, 0, 0, -6, 6, -4, -3, 3, 2, 2, 11, -10, -10, -12, 12, 17, 17, -3, 3, 5, 5, 3, -3, -4, -4, -16]);
ECSearch("36.3-e",(10*w - 56)*OK,[1, -4, -1, 5, -2, 3, 6, -9, -8, 7, 10, -10, 5, -11, -1, 7, -4, 4, 16, -11, 4, -9, 9, 2, -4, 6, 18, 16]);
ECSearch("36.3-f",(10*w - 56)*OK,[-1, -4, -1, -5, 2, -3, -6, 9, -8, 7, -10, 10, 5, 11, 1, -7, 4, 4, 16, -11, 4, -9, 9, 2, -4, -6, -18, -16]);
ECSearch("40.1-a",(-6*w + 34)*OK,[-2, 0, 2, -6, -4, 6, -6, 4, -6, 2, -8, -6, -2, 2, 14, -16, -10, 18, -2, -2, -10, -18, -2, 10, -6, -22, 2, 2]);
ECSearch("40.1-b",(-6*w + 34)*OK,[2, 0, 2, 6, 4, -6, 6, -4, -6, 2, 8, 6, -2, -2, -14, 16, 10, 18, -2, -2, -10, -18, -2, 10, -6, 22, -2, -2]);
ECSearch("40.2-a",(-6*w - 34)*OK,[0, -2, 2, -4, -6, -6, 6, 4, 2, -6, -6, -8, -2, 14, 2, -10, -16, -2, 18, -10, -2, -2, -18, -6, 10, 2, -22, -4]);
ECSearch("40.2-b",(-6*w - 34)*OK,[0, 2, 2, 4, 6, 6, -6, -4, 2, -6, 6, 8, -2, -14, -2, 10, 16, -2, 18, -10, -2, -2, -18, -6, 10, -2, 22, 4]);
ECSearch("44.1-a",(18*w - 100)*OK,[2, 2, 3, 3, 0, 6, -6, 4, 3, -9, -2, -2, -1, 10, 10, 0, 12, -1, -1, 3, -9, -13, -1, 15, -9, 4, -8, 4]);
ECSearch("44.1-b",(18*w - 100)*OK,[2, 0, -3, -3, -2, 8, 0, 2, -5, 5, 8, 10, 9, -6, 12, -14, 4, -9, 3, 17, -17, -11, 7, -15, -15, -10, 12, 10]);
ECSearch("44.1-c",(18*w - 100)*OK,[0, 2, 1, -3, 2, 2, 2, 10, -9, -3, 6, 0, -3, 4, -2, 18, 0, -13, -13, -15, -5, -11, -1, 21, -3, 6, -16, -10]);
ECSearch("44.1-d",(18*w - 100)*OK,[-2, -2, 3, 3, 0, -6, 6, -4, 3, -9, 2, 2, -1, -10, -10, 0, -12, -1, -1, 3, -9, -13, -1, 15, -9, -4, 8, -4]);
ECSearch("44.1-e",(18*w - 100)*OK,[0, -2, 1, -3, -2, -2, -2, -10, -9, -3, -6, 0, -3, -4, 2, -18, 0, -13, -13, -15, -5, -11, -1, 21, -3, -6, 16, 10]);
ECSearch("44.1-f",(18*w - 100)*OK,[-2, 0, -3, -3, 2, -8, 0, -2, -5, 5, -8, -10, 9, 6, -12, 14, -4, -9, 3, 17, -17, -11, 7, -15, -15, 10, -12, -10]);
ECSearch("44.2-a",(-18*w - 100)*OK,[2, 2, 3, 3, 0, -6, 6, 4, -9, 3, -2, -2, -1, 10, 10, 12, 0, -1, -1, -9, 3, -1, -13, -9, 15, -8, 4, 4]);
ECSearch("44.2-b",(-18*w - 100)*OK,[0, 2, -3, -3, -2, 0, 8, 2, 5, -5, 10, 8, 9, 12, -6, 4, -14, 3, -9, -17, 17, 7, -11, -15, -15, 12, -10, -2]);
ECSearch("44.2-c",(-18*w - 100)*OK,[2, 0, -3, 1, 2, 2, 2, 10, -3, -9, 0, 6, -3, -2, 4, 0, 18, -13, -13, -5, -15, -1, -11, -3, 21, -16, 6, -22]);
ECSearch("44.2-d",(-18*w - 100)*OK,[-2, -2, 3, 3, 0, 6, -6, -4, -9, 3, 2, 2, -1, -10, -10, -12, 0, -1, -1, -9, 3, -1, -13, -9, 15, 8, -4, -4]);
ECSearch("44.2-e",(-18*w - 100)*OK,[-2, 0, -3, 1, -2, -2, -2, -10, -3, -9, 0, -6, -3, 2, -4, 0, -18, -13, -13, -5, -15, -1, -11, -3, 21, 16, -6, 22]);
ECSearch("44.2-f",(-18*w - 100)*OK,[0, -2, -3, -3, 2, 0, -8, -2, 5, -5, -10, -8, 9, -12, 6, -4, 14, 3, -9, -17, 17, 7, -11, -15, -15, -12, 10, 2]);
ECSearch("45.1-a",(-3*w - 18)*OK,[1, 4, 0, -4, -2, 0, -8, -2, 8, 4, 6, -2, 16, 8, -2, 0, 12, 10, 18, -6, -6, -10, -4, 18, 8, 22, -6]);
ECSearch("45.1-b",(-3*w - 18)*OK,[-2, -2, 0, -2, 2, 0, -10, 4, 2, -4, 6, -2, 8, 16, -4, 6, -6, -8, 12, -18, -18, 14, 14, -6, -8, 8, 18]);
ECSearch("45.1-c",(-3*w - 18)*OK,[-1, -2, -4, -4, -8, 8, 8, -6, 2, 4, 4, 2, 16, 0, -12, -12, -6, 2, 6, -10, 14, -2, -6, -14, -8, -8, 20]);
ECSearch("45.1-d",(-3*w - 18)*OK,[0, 2, 4, 0, 0, 4, -2, -2, 2, 10, -10, -2, 10, -6, -16, -4, 16, 12, 2, 6, 20, 4, 6, 14, 8, -16, -2]);
ECSearch("45.1-e",(-3*w - 18)*OK,[-2, -2, 0, 2, -2, 0, 10, 4, 2, 4, -6, -2, -8, -16, 4, -6, -6, -8, 12, -18, -18, 14, 14, -6, 8, -8, -18]);
ECSearch("45.1-f",(-3*w - 18)*OK,[-1, -2, 4, 4, 8, -8, -8, -6, 2, -4, -4, 2, -16, 0, 12, 12, -6, 2, 6, -10, 14, -2, -6, -14, 8, 8, -20]);
ECSearch("45.1-g",(-3*w - 18)*OK,[0, 2, -4, 0, 0, -4, 2, -2, 2, -10, 10, -2, -10, 6, 16, 4, 16, 12, 2, 6, 20, 4, 6, 14, -8, 16, 2]);
ECSearch("45.1-h",(-3*w - 18)*OK,[1, 4, 0, 4, 2, 0, 8, -2, 8, -4, -6, -2, -16, -8, 2, 0, 12, 10, 18, -6, -6, -10, -4, 18, -8, -22, 6]);
ECSearch("45.3-a",(19*w - 106)*OK,[0, 3, 1, 1, 6, 8, 6, 3, -8, -6, 4, -3, -6, 2, 2, -4, 18, -16, 5, 12, -19, -12, 2, -1, -6, 0, -7, -11]);
ECSearch("45.3-b",(19*w - 106)*OK,[0, 1, -3, 3, 0, 6, -6, 5, 6, -6, 8, 11, -4, -16, -4, 6, -6, 14, -7, -18, -3, -4, -10, -9, -6, 20, 11, -13]);
ECSearch("45.3-c",(19*w - 106)*OK,[-1, 2, -2, -4, -4, -2, -4, 2, -6, 2, -2, -8, 2, 4, -12, -6, 12, 18, 2, -6, 2, 14, -2, -6, -2, -2, -8, -16]);
ECSearch("45.3-d",(19*w - 106)*OK,[0, -1, -3, -3, 0, -6, 6, -5, 6, -6, -8, -11, -4, 16, 4, -6, 6, 14, -7, -18, -3, -4, -10, -9, -6, -20, -11, 13]);
ECSearch("45.3-e",(19*w - 106)*OK,[0, -3, 1, -1, -6, -8, -6, -3, -8, -6, -4, 3, -6, -2, -2, 4, -18, -16, 5, 12, -19, -12, 2, -1, -6, 0, 7, 11]);
ECSearch("45.3-f",(19*w - 106)*OK,[-1, -2, -2, 4, 4, 2, 4, -2, -6, 2, 2, 8, 2, -4, 12, 6, -12, 18, 2, -6, 2, 14, -2, -6, -2, 2, 8, 16]);
ECSearch("45.4-a",(-2*w + 13)*OK,[-1, 1, 1, 5, -4, 1, -4, -7, 0, 5, 4, 4, -13, 1, 3, 9, -6, 0, 14, 9, -4, 5, 7, -6, -14, 16, 4, -16]);
ECSearch("45.4-b",(-2*w + 13)*OK,[-2, 2, 2, 2, -4, 4, -4, 4, 0, 10, -4, -4, -10, 8, 0, 0, 0, -6, -16, 0, -2, 14, -2, -18, -10, -22, -16, -8]);
ECSearch("45.4-c",(-2*w + 13)*OK,[1, 2, -4, 2, -4, -8, -4, -8, -6, 10, -4, -4, -10, 8, 0, -12, 12, 12, 2, -6, -2, 8, 16, -18, 8, 14, -16, -20]);
ECSearch("45.4-d",(-2*w + 13)*OK,[1, 1, -1, 1, -2, -1, 4, 11, 6, 7, -2, -8, 5, -11, 15, -9, -6, 12, 14, -15, 4, 11, -11, 6, -4, -8, 10, 2]);
ECSearch("45.4-e",(-2*w + 13)*OK,[-2, -2, 2, -2, 4, -4, 4, -4, 0, 10, 4, 4, -10, -8, 0, 0, 0, -6, -16, 0, -2, 14, -2, -18, -10, 22, 16, 8]);
ECSearch("45.4-f",(-2*w + 13)*OK,[1, -1, -1, -1, 2, 1, -4, -11, 6, 7, 2, 8, 5, 11, -15, 9, 6, 12, 14, -15, 4, 11, -11, 6, -4, 8, -10, -2]);
ECSearch("45.4-g",(-2*w + 13)*OK,[1, -2, -4, -2, 4, 8, 4, 8, -6, 10, 4, 4, -10, -8, 0, 12, -12, 12, 2, -6, -2, 8, 16, -18, 8, -14, 16, 20]);
ECSearch("45.4-h",(-2*w + 13)*OK,[-1, -1, 1, -5, 4, -1, 4, 7, 0, 5, -4, -4, -13, -1, -3, -9, 6, 0, 14, 9, -4, 5, 7, -6, -14, -16, -4, 16]);
ECSearch("45.2-a",(3*w - 18)*OK,[1, 4, -4, 0, 0, -2, -8, 8, -2, 6, 4, -2, 8, 16, 0, -2, 10, 12, -6, 18, -10, -6, 18, -4, 22, 8, -10]);
ECSearch("45.2-b",(3*w - 18)*OK,[-2, -2, -2, 0, 0, 2, -10, 2, 4, 6, -4, -2, 16, 8, 6, -4, -8, -6, -18, 12, 14, -18, -6, 14, 8, -8, -14]);
ECSearch("45.2-c",(3*w - 18)*OK,[-1, -2, -4, -4, 8, -8, 8, 2, -6, 4, 4, 2, 0, 16, -12, -12, 2, -6, -10, 6, -2, 14, -14, -6, -8, -8, -12]);
ECSearch("45.2-d",(3*w - 18)*OK,[0, 2, 0, 4, 4, 0, -2, 2, -2, -10, 10, -2, -6, 10, -4, -16, 12, 16, 6, 2, 4, 20, 14, 6, -16, 8, -2]);
ECSearch("45.2-e",(3*w - 18)*OK,[-2, -2, 2, 0, 0, -2, 10, 2, 4, -6, 4, -2, -16, -8, -6, 4, -8, -6, -18, 12, 14, -18, -6, 14, -8, 8, 14]);
ECSearch("45.2-f",(3*w - 18)*OK,[-1, -2, 4, 4, -8, 8, -8, 2, -6, -4, -4, 2, 0, -16, 12, 12, 2, -6, -10, 6, -2, 14, -14, -6, 8, 8, 12]);
ECSearch("45.2-g",(3*w - 18)*OK,[0, 2, 0, -4, -4, 0, 2, 2, -2, 10, -10, -2, 6, -10, 4, 16, 12, 16, 6, 2, 4, 20, 14, 6, 16, -8, 2]);
ECSearch("45.2-h",(3*w - 18)*OK,[1, 4, 4, 0, 0, 2, 8, 8, -2, -6, -4, -2, -8, -16, 0, 2, 10, 12, -6, 18, -10, -6, 18, -4, -22, -8, 10]);
ECSearch("45.6-a",(-19*w - 106)*OK,[0, 3, 1, 6, 1, 6, 8, 3, -6, -8, -3, 4, -6, 2, 2, 18, -4, 5, -16, -19, 12, 2, -12, -6, -1, -7, 0, -4]);
ECSearch("45.6-b",(-19*w - 106)*OK,[0, 1, -3, 0, 3, -6, 6, 5, -6, 6, 11, 8, -4, -4, -16, -6, 6, -7, 14, -3, -18, -10, -4, -6, -9, 11, 20, -4]);
ECSearch("45.6-c",(-19*w - 106)*OK,[-1, 2, -2, -4, -4, -4, -2, 2, 2, -6, -8, -2, 2, -12, 4, 12, -6, 2, 18, 2, -6, -2, 14, -2, -6, -8, -2, 12]);
ECSearch("45.6-d",(-19*w - 106)*OK,[0, -1, -3, 0, -3, 6, -6, -5, -6, 6, -11, -8, -4, 4, 16, 6, -6, -7, 14, -3, -18, -10, -4, -6, -9, -11, -20, 4]);
ECSearch("45.6-e",(-19*w - 106)*OK,[0, -3, 1, -6, -1, -6, -8, -3, -6, -8, 3, -4, -6, -2, -2, -18, 4, 5, -16, -19, 12, 2, -12, -6, -1, 7, 0, 4]);
ECSearch("45.6-f",(-19*w - 106)*OK,[-1, -2, -2, 4, 4, 4, 2, -2, 2, -6, 8, 2, 2, 12, -4, -12, 6, 2, 18, 2, -6, -2, 14, -2, -6, 8, 2, -12]);
ECSearch("45.5-a",(2*w + 13)*OK,[-1, 1, 1, -4, 5, -4, 1, -7, 5, 0, 4, 4, -13, 3, 1, -6, 9, 14, 0, -4, 9, 7, 5, -14, -6, 4, 16, 6]);
ECSearch("45.5-b",(2*w + 13)*OK,[-2, 2, 2, -4, 2, -4, 4, 4, 10, 0, -4, -4, -10, 0, 8, 0, 0, -16, -6, -2, 0, -2, 14, -10, -18, -16, -22, 12]);
ECSearch("45.5-c",(2*w + 13)*OK,[1, 2, -4, -4, 2, -4, -8, -8, 10, -6, -4, -4, -10, 0, 8, 12, -12, 2, 12, -2, -6, 16, 8, 8, -18, -16, 14, -12]);
ECSearch("45.5-d",(2*w + 13)*OK,[1, 1, -1, -2, 1, 4, -1, 11, 7, 6, -8, -2, 5, 15, -11, -6, -9, 14, 12, 4, -15, -11, 11, -4, 6, 10, -8, 0]);
ECSearch("45.5-e",(2*w + 13)*OK,[-2, -2, 2, 4, -2, 4, -4, -4, 10, 0, 4, 4, -10, 0, -8, 0, 0, -16, -6, -2, 0, -2, 14, -10, -18, 16, 22, -12]);
ECSearch("45.5-f",(2*w + 13)*OK,[1, -1, -1, 2, -1, -4, 1, -11, 7, 6, 8, 2, 5, -15, 11, 6, 9, 14, 12, 4, -15, -11, 11, -4, 6, -10, 8, 0]);
ECSearch("45.5-g",(2*w + 13)*OK,[1, -2, -4, 4, -2, 4, 8, 8, 10, -6, 4, 4, -10, 0, -8, -12, 12, 2, 12, -2, -6, 16, 8, 8, -18, 16, -14, 12]);
ECSearch("45.5-h",(2*w + 13)*OK,[-1, -1, 1, 4, -5, 4, -1, 7, 5, 0, -4, -4, -13, -3, -1, 6, -9, 14, 0, -4, 9, 7, 5, -14, -6, -4, -16, -6]);
ECSearch("48.1-a",(8*w + 44)*OK,[1, -2, 3, -3, 2, 7, -8, -1, -10, -5, -2, 8, -3, 5, 15, -9, 6, 16, -14, 17, -8, 15, 5, 6, 16, 8, 18, 22]);
ECSearch("48.1-b",(8*w + 44)*OK,[0, 3, 3, -4, -4, 4, 0, -4, -5, 3, -4, 8, -9, 0, 4, 8, -16, -13, 3, 7, 7, -5, 19, -5, 3, -12, -4, -8]);
ECSearch("48.1-c",(8*w + 44)*OK,[1, -2, -3, -3, 2, 1, 4, -7, 2, -5, -2, -4, -3, -13, 9, -9, 6, 4, -2, -1, -8, 9, 11, -18, -8, -4, -18, -14]);
ECSearch("48.1-d",(8*w + 44)*OK,[0, 3, -1, 0, 4, 4, -4, -8, -1, 3, 8, -8, -1, -12, -8, -12, -12, 7, -17, -13, 15, -17, 15, 19, -1, 16, 4, 0]);
ECSearch("48.1-e",(8*w + 44)*OK,[3, 0, 3, 1, -2, -1, 0, 1, 10, -3, 4, 10, -3, -3, 11, -11, -2, 8, -12, 7, -2, 7, 7, -14, -6, -18, -8, 20]);
ECSearch("48.1-f",(8*w + 44)*OK,[-3, 0, 3, -1, 2, 1, 0, -1, 10, -3, -4, -10, -3, 3, -11, 11, 2, 8, -12, 7, -2, 7, 7, -14, -6, 18, 8, -20]);
ECSearch("48.1-g",(8*w + 44)*OK,[0, 3, -1, 0, -4, -4, 4, 8, -1, 3, -8, 8, -1, 12, 8, 12, 12, 7, -17, -13, 15, -17, 15, 19, -1, -16, -4, 0]);
ECSearch("48.1-h",(8*w + 44)*OK,[-1, -2, -3, 3, -2, -1, -4, 7, 2, -5, 2, 4, -3, 13, -9, 9, -6, 4, -2, -1, -8, 9, 11, -18, -8, 4, 18, 14]);
ECSearch("48.1-i",(8*w + 44)*OK,[0, 3, 3, 4, 4, -4, 0, 4, -5, 3, 4, -8, -9, 0, -4, -8, 16, -13, 3, 7, 7, -5, 19, -5, 3, 12, 4, 8]);
ECSearch("48.1-j",(8*w + 44)*OK,[-1, -2, 3, 3, -2, -7, 8, 1, -10, -5, 2, -8, -3, -5, -15, 9, -6, 16, -14, 17, -8, 15, 5, 6, 16, -8, -18, -22]);
ECSearch("48.2-a",(8*w - 44)*OK,[1, 3, -2, 2, -3, -8, 7, -1, -5, -10, 8, -2, -3, 15, 5, 6, -9, -14, 16, -8, 17, 5, 15, 16, 6, 18, 8, -8]);
ECSearch("48.2-b",(8*w - 44)*OK,[0, 3, 3, -4, -4, 0, 4, -4, 3, -5, 8, -4, -9, 4, 0, -16, 8, 3, -13, 7, 7, 19, -5, 3, -5, -4, -12, -16]);
ECSearch("48.2-c",(8*w - 44)*OK,[1, -3, -2, 2, -3, 4, 1, -7, -5, 2, -4, -2, -3, 9, -13, 6, -9, -2, 4, -8, -1, 11, 9, -8, -18, -18, -4, 4]);
ECSearch("48.2-d",(8*w - 44)*OK,[0, -1, 3, 4, 0, -4, 4, -8, 3, -1, -8, 8, -1, -8, -12, -12, -12, -17, 7, 15, -13, 15, -17, -1, 19, 4, 16, 4]);
ECSearch("48.2-e",(8*w - 44)*OK,[3, 3, 0, -2, 1, 0, -1, 1, -3, 10, 10, 4, -3, 11, -3, -2, -11, -12, 8, -2, 7, 7, 7, -6, -14, -8, -18, 10]);
ECSearch("48.2-f",(8*w - 44)*OK,[-3, 3, 0, 2, -1, 0, 1, -1, -3, 10, -10, -4, -3, -11, 3, 2, 11, -12, 8, -2, 7, 7, 7, -6, -14, 8, 18, -10]);
ECSearch("48.2-g",(8*w - 44)*OK,[0, -1, 3, -4, 0, 4, -4, 8, 3, -1, 8, -8, -1, 8, 12, 12, 12, -17, 7, 15, -13, 15, -17, -1, 19, -4, -16, -4]);
ECSearch("48.2-h",(8*w - 44)*OK,[-1, -3, -2, -2, 3, -4, -1, 7, -5, 2, 4, 2, -3, -9, 13, -6, 9, -2, 4, -8, -1, 11, 9, -8, -18, 18, 4, -4]);
ECSearch("48.2-i",(8*w - 44)*OK,[0, 3, 3, 4, 4, 0, -4, 4, 3, -5, -8, 4, -9, -4, 0, 16, -8, 3, -13, 7, 7, 19, -5, 3, -5, 4, 12, 16]);
ECSearch("48.2-j",(8*w - 44)*OK,[-1, 3, -2, -2, 3, 8, -7, 1, -5, -10, -8, 2, -3, -15, -5, -6, 9, -14, 16, -8, 17, 5, 15, 16, 6, -18, -8, 8]);
ECSearch("49.1-a",(7)*OK,[1, 0, 0, -2, -2, -6, 6, -6, 6, 0, 10, 10, 6, -6, -6, 6, -12, 12, 10, 10, -10, -10, 2, 2, -10, -10, -6, 6, 0]);
ECSearch("49.1-b",(7)*OK,[1, 0, 0, -2, -2, 6, -6, 6, -6, 0, 10, 10, -6, 6, 6, -6, 12, -12, 10, 10, -10, -10, 2, 2, -10, -10, 6, -6, 0]);
ECSearch("50.2-a",(w + 9)*OK,[-1, 1, 0, 0, -3, 6, 0, 2, 6, 3, -1, 7, 1, 10, -14, -12, -15, 2, 10, -6, -12, -2, -2, -21, -3, 14, 4, -4]);
ECSearch("50.2-b",(w + 9)*OK,[0, 1, 1, -6, -3, -7, 0, -8, -4, 10, -5, 0, -6, -1, 11, 4, -4, 9, -15, 3, 0, 20, 1, 8, 12, -15, 0, 14]);
ECSearch("50.2-c",(w + 9)*OK,[1, -1, 0, 0, 3, -6, 0, -2, 6, 3, 1, -7, 1, -10, 14, 12, 15, 2, 10, -6, -12, -2, -2, -21, -3, -14, -4, 4]);
ECSearch("50.2-d",(w + 9)*OK,[0, -1, 1, 6, 3, 7, 0, 8, -4, 10, 5, 0, -6, 1, -11, -4, 4, 9, -15, 3, 0, 20, 1, 8, 12, 15, 0, -14]);
ECSearch("50.3-a",(w - 9)*OK,[1, -1, 0, -3, 0, 0, 6, 2, 3, 6, 7, -1, 1, -14, 10, -15, -12, 10, 2, -12, -6, -2, -2, -3, -21, 4, 14, -19]);
ECSearch("50.3-b",(w - 9)*OK,[1, 0, 1, -3, -6, 0, -7, -8, 10, -4, 0, -5, -6, 11, -1, -4, 4, -15, 9, 0, 3, 1, 20, 12, 8, 0, -15, -21]);
ECSearch("50.3-c",(w - 9)*OK,[-1, 1, 0, 3, 0, 0, -6, -2, 3, 6, -7, 1, 1, 14, -10, 15, 12, 10, 2, -12, -6, -2, -2, -3, -21, -4, -14, 19]);
ECSearch("50.3-d",(w - 9)*OK,[-1, 0, 1, 3, 6, 0, 7, 8, 10, -4, 0, 5, -6, -11, 1, 4, -4, -15, 9, 0, 3, 1, 20, 12, 8, 0, 15, 21]);
ECSearch("54.1-a",(-3*w + 15)*OK,[1, 2, -5, 0, -1, 0, -7, -3, 5, -9, 2, -1, 4, 14, 4, -9, -13, 2, 3, 6, 11, -4, -10, -15, -14, 12, -10]);
ECSearch("54.1-b",(-3*w + 15)*OK,[3, 3, -3, 3, -6, -6, -8, -12, 0, -8, 4, 2, -17, 1, -9, 9, 11, -13, 6, -18, -4, -16, 12, -12, 16, -8, 4]);
ECSearch("54.1-c",(-3*w + 15)*OK,[0, -3, -6, -3, 0, -3, 7, -3, 9, -2, 1, -13, 10, 10, -3, -6, 14, -1, 0, 15, -4, 11, 21, 6, -2, -8, 4]);
ECSearch("54.1-d",(-3*w + 15)*OK,[3, 0, -3, -6, -3, 0, -7, -9, 3, -1, 2, -13, -10, -10, -6, -3, -1, 14, -15, 0, 11, -4, -6, -21, 8, 2, 20]);
ECSearch("54.1-e",(-3*w + 15)*OK,[1, 0, -3, 2, 1, 8, -1, 1, 1, 1, -10, 3, -10, -6, -18, 9, -13, 14, 7, 16, 15, 4, -18, 5, -8, 18, 4]);
ECSearch("54.1-f",(-3*w + 15)*OK,[-2, -1, -4, 3, 4, -3, 1, 3, 5, 6, -5, 11, -10, 10, 17, 6, -10, -7, -18, 15, -16, 17, 5, -6, -10, -6, 4]);
ECSearch("54.1-g",(-3*w + 15)*OK,[-3, 2, -2, -3, 4, -1, -7, 5, 2, -4, -2, -3, 9, -13, -6, -9, 2, -4, 8, -1, -11, -9, -8, 18, 18, 4, 4]);
ECSearch("54.1-h",(-3*w + 15)*OK,[-1, 2, -1, 4, -9, 0, -9, 7, -5, 1, -10, -1, 4, 2, 4, -1, 7, 10, 1, -2, 15, -12, 2, -13, 18, -12, -2]);
ECSearch("54.1-i",(-3*w + 15)*OK,[0, -1, -2, 3, -8, -1, -1, -1, -1, -10, 1, 3, -6, -10, -9, 18, 14, -13, -16, -7, 4, 15, -5, 18, 18, -8, -20]);
ECSearch("54.1-j",(-3*w + 15)*OK,[3, -1, -1, -3, -4, -2, 4, 8, -10, 4, -10, 6, -15, -5, 3, 9, 5, -7, 2, -10, 4, -18, 10, -6, 0, -4, -4]);
ECSearch("54.1-k",(-3*w + 15)*OK,[-3, -3, -3, 3, 6, 6, -8, 0, 12, 4, -8, 2, 1, -17, -9, 9, -13, 11, 18, -6, -16, -4, 12, -12, -8, 16, 16]);
ECSearch("54.1-l",(-3*w + 15)*OK,[-1, 2, 1, -4, 9, 0, 9, 7, -5, -1, 10, -1, -4, -2, -4, 1, 7, 10, 1, -2, 15, -12, 2, -13, -18, 12, 2]);
ECSearch("54.1-m",(-3*w + 15)*OK,[3, 3, 3, -3, 6, 6, 8, -12, 0, 8, -4, 2, 17, -1, 9, -9, 11, -13, 6, -18, -4, -16, 12, -12, -16, 8, -4]);
ECSearch("54.1-n",(-3*w + 15)*OK,[3, -1, 1, 3, 4, 2, -4, 8, -10, -4, 10, 6, 15, 5, -3, -9, 5, -7, 2, -10, 4, -18, 10, -6, 0, 4, 4]);
ECSearch("54.1-o",(-3*w + 15)*OK,[-3, 2, 2, 3, -4, 1, 7, 5, 2, 4, 2, -3, -9, 13, 6, 9, 2, -4, 8, -1, -11, -9, -8, 18, -18, -4, -4]);
ECSearch("54.1-p",(-3*w + 15)*OK,[3, 0, 3, 6, 3, 0, 7, -9, 3, 1, -2, -13, 10, 10, 6, 3, -1, 14, -15, 0, 11, -4, -6, -21, -8, -2, -20]);
ECSearch("54.1-q",(-3*w + 15)*OK,[0, -1, 2, -3, 8, 1, 1, -1, -1, 10, -1, 3, 6, 10, 9, -18, 14, -13, -16, -7, 4, 15, -5, 18, -18, 8, 20]);
ECSearch("54.1-r",(-3*w + 15)*OK,[1, 0, 3, -2, -1, -8, 1, 1, 1, -1, 10, 3, 10, 6, 18, -9, -13, 14, 7, 16, 15, 4, -18, 5, 8, -18, -4]);
ECSearch("54.1-s",(-3*w + 15)*OK,[1, 2, 5, 0, 1, 0, 7, -3, 5, 9, -2, -1, -4, -14, -4, 9, -13, 2, 3, 6, 11, -4, -10, -15, 14, -12, 10]);
ECSearch("54.1-t",(-3*w + 15)*OK,[-2, -1, 4, -3, -4, 3, -1, 3, 5, -6, 5, 11, 10, -10, -17, -6, -10, -7, -18, 15, -16, 17, 5, -6, 10, 6, -4]);
ECSearch("54.1-u",(-3*w + 15)*OK,[-3, -3, 3, -3, -6, -6, 8, 0, 12, -4, 8, 2, -1, 17, 9, -9, -13, 11, 18, -6, -16, -4, 12, -12, 8, -16, -16]);
ECSearch("54.1-v",(-3*w + 15)*OK,[0, -3, 6, 3, 0, 3, -7, -3, 9, 2, -1, -13, -10, -10, 3, 6, 14, -1, 0, 15, -4, 11, 21, 6, 2, 8, -4]);
ECSearch("54.2-a",(3*w + 15)*OK,[2, 1, 0, -5, 0, -1, -7, 5, -3, 2, -9, -1, 14, 4, -9, 4, 2, -13, 6, 3, -4, 11, -15, -10, 12, -14, 12]);
ECSearch("54.2-b",(3*w + 15)*OK,[3, 3, 3, -3, -6, -6, -8, 0, -12, 4, -8, 2, 1, -17, 9, -9, -13, 11, -18, 6, -16, -4, -12, 12, -8, 16, 16]);
ECSearch("54.2-c",(3*w + 15)*OK,[-3, 0, -3, -6, -3, 0, 7, 9, -3, 1, -2, -13, 10, 10, -6, -3, -1, 14, 15, 0, 11, -4, 6, 21, -8, -2, -20]);
ECSearch("54.2-d",(3*w + 15)*OK,[0, 3, -6, -3, 0, -3, -7, 3, -9, 2, -1, -13, -10, -10, -3, -6, 14, -1, 0, -15, -4, 11, -21, -6, 2, 8, -4]);
ECSearch("54.2-e",(3*w + 15)*OK,[0, 1, 2, -3, 8, 1, -1, 1, 1, -10, 1, 3, -6, -10, 9, -18, 14, -13, 16, 7, 4, 15, 5, -18, 18, -8, -20]);
ECSearch("54.2-f",(3*w + 15)*OK,[-1, -2, 3, -4, -3, 4, 1, 5, 3, -5, 6, 11, 10, -10, 6, 17, -7, -10, 15, -18, 17, -16, -6, 5, -6, -10, 6]);
ECSearch("54.2-g",(3*w + 15)*OK,[2, -3, -3, -2, -1, 4, -7, 2, 5, -2, -4, -3, -13, 9, -9, -6, -4, 2, -1, 8, -9, -11, 18, -8, 4, 18, -14]);
ECSearch("54.2-h",(3*w + 15)*OK,[2, -1, 4, -1, 0, -9, -9, -5, 7, -10, 1, -1, 2, 4, -1, 4, 10, 7, -2, 1, -12, 15, -13, 2, -12, 18, -4]);
ECSearch("54.2-i",(3*w + 15)*OK,[-1, 0, 3, -2, -1, -8, -1, -1, -1, 1, -10, 3, -10, -6, 18, -9, -13, 14, -7, -16, 15, 4, 18, -5, -8, 18, 4]);
ECSearch("54.2-j",(3*w + 15)*OK,[-1, 3, -3, -1, -2, -4, 4, -10, 8, -10, 4, 6, -5, -15, 9, 3, -7, 5, -10, 2, -18, 4, -6, 10, -4, 0, 14]);
ECSearch("54.2-k",(3*w + 15)*OK,[-3, -3, 3, -3, 6, 6, -8, 12, 0, -8, 4, 2, -17, 1, 9, -9, 11, -13, -6, 18, -4, -16, -12, 12, 16, -8, 4]);
ECSearch("54.2-l",(3*w + 15)*OK,[2, -1, -4, 1, 0, 9, 9, -5, 7, 10, -1, -1, -2, -4, 1, -4, 10, 7, -2, 1, -12, 15, -13, 2, 12, -18, 4]);
ECSearch("54.2-m",(3*w + 15)*OK,[3, 3, -3, 3, 6, 6, 8, 0, -12, -4, 8, 2, -1, 17, -9, 9, -13, 11, -18, 6, -16, -4, -12, 12, 8, -16, -16]);
ECSearch("54.2-n",(3*w + 15)*OK,[-1, 3, 3, 1, 2, 4, -4, -10, 8, 10, -4, 6, 5, 15, -9, -3, -7, 5, -10, 2, -18, 4, -6, 10, 4, 0, -14]);
ECSearch("54.2-o",(3*w + 15)*OK,[2, -3, 3, 2, 1, -4, 7, 2, 5, 2, 4, -3, 13, -9, 9, 6, -4, 2, -1, 8, -9, -11, 18, -8, -4, -18, 14]);
ECSearch("54.2-p",(3*w + 15)*OK,[0, 3, 6, 3, 0, 3, 7, 3, -9, -2, 1, -13, 10, 10, 3, 6, 14, -1, 0, -15, -4, 11, -21, -6, -2, -8, 4]);
ECSearch("54.2-q",(3*w + 15)*OK,[-1, 0, -3, 2, 1, 8, 1, -1, -1, -1, 10, 3, 10, 6, -18, 9, -13, 14, -7, -16, 15, 4, 18, -5, 8, -18, -4]);
ECSearch("54.2-r",(3*w + 15)*OK,[0, 1, -2, 3, -8, -1, 1, 1, 1, 10, -1, 3, 6, 10, -9, 18, 14, -13, 16, 7, 4, 15, 5, -18, -18, 8, 20]);
ECSearch("54.2-s",(3*w + 15)*OK,[2, 1, 0, 5, 0, 1, 7, 5, -3, -2, 9, -1, -14, -4, 9, -4, 2, -13, 6, 3, -4, 11, -15, -10, -12, 14, -12]);
ECSearch("54.2-t",(3*w + 15)*OK,[-1, -2, -3, 4, 3, -4, -1, 5, 3, 5, -6, 11, -10, 10, -6, -17, -7, -10, 15, -18, 17, -16, -6, 5, 6, 10, -6]);
ECSearch("54.2-u",(3*w + 15)*OK,[-3, -3, -3, 3, -6, -6, 8, 12, 0, 8, -4, 2, 17, -1, -9, 9, 11, -13, -6, 18, -4, -16, -12, 12, -16, 8, -4]);
ECSearch("54.2-v",(3*w + 15)*OK,[-3, 0, 3, 6, 3, 0, -7, 9, -3, -1, 2, -13, -10, -10, 6, 3, -1, 14, 15, 0, 11, -4, 6, 21, 8, 2, 20]);
ECSearch("55.1-a",(13*w + 72)*OK,[-1, 3, 2, 3, 2, 1, 2, -4, 0, -9, 0, 8, 2, 5, -4, -9, -1, 0, 2, 0, -3, -17, 10, 8, 11, -14, 4, -16]);
ECSearch("55.1-b",(13*w + 72)*OK,[-1, -3, -2, 3, -2, -1, -2, 4, 0, -9, 0, -8, 2, -5, 4, 9, 1, 0, 2, 0, -3, -17, 10, 8, 11, 14, -4, 16]);
ECSearch("55.4-a",(-13*w + 72)*OK,[-1, 2, 3, 3, 2, 2, 1, -4, -9, 0, 8, 0, 2, -4, 5, -1, -9, 2, 0, -3, 0, 10, -17, 11, 8, 4, -14, -13]);
ECSearch("55.4-b",(-13*w + 72)*OK,[-1, -2, -3, 3, -2, -2, -1, 4, -9, 0, -8, 0, 2, 4, -5, 1, 9, 2, 0, -3, 0, 10, -17, 11, 8, -4, 14, 13]);