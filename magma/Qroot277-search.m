print "Field 2.2.277.1";
Qx<x> := PolynomialRing(RationalField());
K<w> := NumberField(x^2 - x - 69);
OK := Integers(K);
Plist := [];
Append(~Plist,(-w + 9)*OK);
Append(~Plist,(w + 8)*OK);
Append(~Plist,(2)*OK);
Append(~Plist,(-6*w - 47)*OK);
Append(~Plist,(6*w - 53)*OK);
Append(~Plist,(w - 8)*OK);
Append(~Plist,(w + 7)*OK);
Append(~Plist,(4*w + 31)*OK);
Append(~Plist,(4*w - 35)*OK);
Append(~Plist,(3*w + 23)*OK);
Append(~Plist,(3*w - 26)*OK);
Append(~Plist,(5)*OK);
Append(~Plist,(7*w + 55)*OK);
Append(~Plist,(7*w - 62)*OK);
Append(~Plist,(w + 10)*OK);
Append(~Plist,(w - 11)*OK);
Append(~Plist,(-2*w + 19)*OK);
Append(~Plist,(-2*w - 17)*OK);
Append(~Plist,(-9*w + 79)*OK);
Append(~Plist,(9*w + 70)*OK);
Append(~Plist,(w + 1)*OK);
Append(~Plist,(w - 2)*OK);
Append(~Plist,(-3*w - 22)*OK);
Append(~Plist,(3*w - 25)*OK);
Append(~Plist,(-3*w + 28)*OK);
Append(~Plist,(-3*w - 25)*OK);
Append(~Plist,(29*w + 227)*OK);
Append(~Plist,(-27*w - 211)*OK);
Append(~Plist,(21*w + 164)*OK);
Append(~Plist,(-21*w + 185)*OK);
Append(~Plist,(w + 13)*OK);
Append(~Plist,(w - 14)*OK);
Append(~Plist,(11)*OK);
Append(~Plist,(10*w + 79)*OK);
Append(~Plist,(10*w - 89)*OK);
Append(~Plist,(26*w + 203)*OK);
Append(~Plist,(-26*w + 229)*OK);
Append(~Plist,(25*w - 221)*OK);
Append(~Plist,(25*w + 196)*OK);
Append(~Plist,(-19*w + 167)*OK);
Append(~Plist,(19*w + 148)*OK);
Append(~Plist,(-6*w + 55)*OK);
Append(~Plist,(-6*w - 49)*OK);
Append(~Plist,(12*w + 95)*OK);
Append(~Plist,(12*w - 107)*OK);
Append(~Plist,(-2*w + 7)*OK);
Append(~Plist,(2*w + 5)*OK);
Append(~Plist,(-2*w + 1)*OK);
Append(~Plist,(-3*w + 20)*OK);
Append(~Plist,(3*w + 17)*OK);
Append(~Plist,(17)*OK);
Append(~Plist,(-17*w + 151)*OK);
Append(~Plist,(17*w + 134)*OK);
Append(~Plist,(-3*w - 29)*OK);
Append(~Plist,(3*w - 32)*OK);
Append(~Plist,(w + 19)*OK);
Append(~Plist,(w - 20)*OK);
Append(~Plist,(-7*w + 59)*OK);
Append(~Plist,(7*w + 52)*OK);
Append(~Plist,(-3*w - 16)*OK);
Append(~Plist,(3*w - 19)*OK);
Append(~Plist,(36*w + 281)*OK);
Append(~Plist,(-36*w + 317)*OK);
Append(~Plist,(-9*w - 68)*OK);
Append(~Plist,(9*w - 77)*OK);
Append(~Plist,(-27*w - 212)*OK);
Append(~Plist,(27*w - 239)*OK);
Append(~Plist,(-4*w - 25)*OK);
Append(~Plist,(4*w - 29)*OK);
Append(~Plist,(-3*w + 17)*OK);
Append(~Plist,(3*w + 14)*OK);
Append(~Plist,(7*w + 58)*OK);
Append(~Plist,(-7*w + 65)*OK);
Append(~Plist,(-9*w + 82)*OK);
Append(~Plist,(-9*w - 73)*OK);
Append(~Plist,(15*w + 119)*OK);
Append(~Plist,(15*w - 134)*OK);
Append(~Plist,(5*w + 44)*OK);
Append(~Plist,(-5*w + 49)*OK);
Append(~Plist,(-3*w + 34)*OK);
Append(~Plist,(3*w + 31)*OK);
Append(~Plist,(-8*w - 59)*OK);
Append(~Plist,(-8*w + 67)*OK);
Append(~Plist,(-3*w + 14)*OK);
Append(~Plist,(3*w + 11)*OK);
Append(~Plist,(43*w - 380)*OK);
Append(~Plist,(43*w + 337)*OK);
Append(~Plist,(-33*w + 292)*OK);
Append(~Plist,(33*w + 259)*OK);
Append(~Plist,(-3*w - 10)*OK);
Append(~Plist,(3*w - 13)*OK);
Append(~Plist,(-3*w - 32)*OK);
Append(~Plist,(3*w - 35)*OK);
Append(~Plist,(-5*w + 37)*OK);
Append(~Plist,(-5*w - 32)*OK);
Append(~Plist,(-11*w - 83)*OK);
Append(~Plist,(11*w - 94)*OK);
Append(~Plist,(-6*w + 47)*OK);
Append(~Plist,(6*w + 41)*OK);
Append(~Plist,(12*w - 103)*OK);
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

ECSearch("3.1-a",(-w + 9)*OK,[-2, -1, -2, 0, -3, 1, 2, -2, 8, 0, 2, -5, 9, -3, -9, 8, 4, -4, 6, -16, 4, -10, -8, -8, -2, 10, 4, 17, 15]);
ECSearch("3.2-a",(w + 8)*OK,[-2, -1, 0, -2, 1, -3, -2, 2, 0, 8, 2, 9, -5, -9, -3, 4, 8, 6, -4, 4, -16, -8, -10, -2, -8, 4, 10, 15, 17]);
ECSearch("4.1-a",(2)*OK,[2, -3, -3, 2, 2, 2, 8, 3, 9, 4, -1, 6, 1, 6, -4, 0, 10, 5, -10, 5, 10, -8, -3, -12, -7, -4, 16, -6, -1]);
ECSearch("4.1-b",(2)*OK,[-3, 2, 2, -3, 2, 2, 3, 8, 4, 9, -1, 1, 6, -4, 6, 10, 0, -10, 5, 10, 5, -3, -8, -7, -12, 16, -4, -1, -6]);
ECSearch("9.1-a",(3)*OK,[-1, 4, -3, -6, 1, 2, -5, -7, 0, -7, -5, 9, 9, 9, 11, 4, 2, 9, -4, 10, 2, -5, 13, -8, 4, 4, 14, 0]);
ECSearch("9.1-b",(3)*OK,[-1, -3, 4, 1, -6, -5, 2, 0, -7, -7, 9, -5, 9, 9, 4, 11, 9, 2, 10, -4, -5, 2, -8, 13, 4, 4, 0, 14]);
ECSearch("9.2-a",(5*w - 44)*OK,[2, 1, 0, 2, 1, -3, -2, 2, 0, -8, -2, -9, -5, 9, -3, -4, 8, 6, 4, 4, -16, 8, -10, -2, -8, 4, -10, 15, -17]);
ECSearch("9.3-a",(5*w + 39)*OK,[2, 1, 2, 0, -3, 1, 2, -2, -8, 0, -2, -5, -9, -3, 9, 8, -4, 4, 6, -16, 4, -10, 8, -8, -2, -10, 4, -17, 15]);
ECSearch("12.1-a",(-2*w + 18)*OK,[-1, 1, 2, 4, 0, 8, -5, -7, 4, -3, 2, 5, 8, 6, 0, 6, -9, 6, -3, 8, -6, 5, 4, 3, 0, -2, -14, 7]);
ECSearch("12.1-b",(-2*w + 18)*OK,[3, 5, -2, 4, 4, -4, 3, -3, 4, 1, -6, 1, 8, -6, -8, 6, -5, 2, -3, 4, 6, -15, 4, 11, -8, 6, -6, 15]);
ECSearch("12.1-c",(-2*w + 18)*OK,[2, -2, 1, -4, 2, 2, -7, 6, -6, 4, 9, 0, 0, 6, 0, 3, 0, 6, -7, 8, 6, 6, 2, -1, -6, 6, 6, -3]);
ECSearch("12.1-d",(-2*w + 18)*OK,[0, 2, -2, -5, -5, 2, -6, 6, -2, -2, -3, 1, -7, -3, -8, 0, -2, 2, -6, -14, -6, 6, 4, -16, -2, -6, -9, 3]);
ECSearch("12.1-e",(-2*w + 18)*OK,[-3, 3, -4, -4, 2, 2, 3, -9, 4, -1, -6, -5, 0, -4, 0, -2, 5, -4, -7, -2, 16, -9, 12, -1, -16, -14, 6, -13]);
ECSearch("12.1-f",(-2*w + 18)*OK,[2, -2, -4, 1, -3, 2, -2, -4, 4, -6, -1, 5, 5, -9, 0, -12, 0, 6, -12, 8, 6, -4, -8, -6, -6, 16, 1, 7]);
ECSearch("12.2-a",(2*w + 16)*OK,[-1, 2, 1, 0, 4, -5, 8, 4, -7, -3, 5, 2, 6, 8, 6, 0, 6, -9, 8, -3, 5, -6, 3, 4, -2, 0, 7, -14]);
ECSearch("12.2-b",(2*w + 16)*OK,[3, -2, 5, 4, 4, 3, -4, 4, -3, 1, 1, -6, -6, 8, 6, -8, 2, -5, 4, -3, -15, 6, 11, 4, 6, -8, 15, -6]);
ECSearch("12.2-c",(2*w + 16)*OK,[2, 1, -2, 2, -4, -7, 2, -6, 6, 4, 0, 9, 6, 0, 3, 0, 6, 0, 8, -7, 6, 6, -1, 2, 6, -6, -3, 6]);
ECSearch("12.2-d",(2*w + 16)*OK,[0, -2, 2, -5, -5, -6, 2, -2, 6, -2, 1, -3, -3, -7, 0, -8, 2, -2, -14, -6, 6, -6, -16, 4, -6, -2, 3, -9]);
ECSearch("12.2-e",(2*w + 16)*OK,[-3, -4, 3, 2, -4, 3, 2, 4, -9, -1, -5, -6, -4, 0, -2, 0, -4, 5, -2, -7, -9, 16, -1, 12, -14, -16, -13, 6]);
ECSearch("12.2-f",(2*w + 16)*OK,[2, -4, -2, -3, 1, -2, 2, 4, -4, -6, 5, -1, -9, 5, -12, 0, 6, 0, 8, -12, -4, 6, -6, -8, 16, -6, 7, 1]);
ECSearch("13.1-a",(w - 8)*OK,[-2, 2, -1, 4, 1, -2, 3, 8, 3, 0, -6, 4, -8, -2, 3, -8, -12, -7, -14, -6, -4, 3, 8, -4, 0, 8, 9, -14, 7]);
ECSearch("13.2-a",(w + 7)*OK,[2, -2, -1, 1, 4, -2, 8, 3, 0, 3, -6, -8, 4, 3, -2, -12, -8, -14, -7, -4, -6, 8, 3, 0, -4, 9, 8, 7, -14]);