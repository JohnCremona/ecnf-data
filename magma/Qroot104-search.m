print "Field 2.2.104.1";
Qx<x> := PolynomialRing(RationalField());
K<w> := NumberField(x^2 - 26);
OK := Integers(K);
Plist := [];
Append(~Plist,(2)*OK+(w)*OK);
Append(~Plist,(5)*OK+(w + 1)*OK);
Append(~Plist,(5)*OK+(w + 4)*OK);
Append(~Plist,(3)*OK);
Append(~Plist,(11)*OK+(w + 2)*OK);
Append(~Plist,(11)*OK+(w + 9)*OK);
Append(~Plist,(13)*OK+(w)*OK);
Append(~Plist,(w + 3)*OK);
Append(~Plist,(w - 3)*OK);
Append(~Plist,(19)*OK+(w + 8)*OK);
Append(~Plist,(19)*OK+(w + 11)*OK);
Append(~Plist,(w + 7)*OK);
Append(~Plist,(w - 7)*OK);
Append(~Plist,(37)*OK+(w + 10)*OK);
Append(~Plist,(37)*OK+(w + 27)*OK);
Append(~Plist,(7)*OK);
Append(~Plist,(59)*OK+(w + 12)*OK);
Append(~Plist,(59)*OK+(w + 47)*OK);
Append(~Plist,(67)*OK+(w + 19)*OK);
Append(~Plist,(67)*OK+(w + 48)*OK);
Append(~Plist,(-2*w + 5)*OK);
Append(~Plist,(2*w + 5)*OK);
Append(~Plist,(83)*OK+(w + 21)*OK);
Append(~Plist,(83)*OK+(w + 62)*OK);
Append(~Plist,(-2*w + 1)*OK);
Append(~Plist,(-2*w - 1)*OK);
Append(~Plist,(109)*OK+(w + 35)*OK);
Append(~Plist,(109)*OK+(w + 74)*OK);
Append(~Plist,(-3*w + 11)*OK);
Append(~Plist,(3*w + 11)*OK);
Append(~Plist,(-3*w + 19)*OK);
Append(~Plist,(3*w + 19)*OK);
Append(~Plist,(149)*OK+(w + 18)*OK);
Append(~Plist,(149)*OK+(w + 131)*OK);
Append(~Plist,(163)*OK+(w + 29)*OK);
Append(~Plist,(163)*OK+(w + 134)*OK);
Append(~Plist,(-4*w + 15)*OK);
Append(~Plist,(4*w + 15)*OK);
Append(~Plist,(197)*OK+(w + 82)*OK);
Append(~Plist,(197)*OK+(w + 115)*OK);
Append(~Plist,(w + 15)*OK);
Append(~Plist,(w - 15)*OK);
Append(~Plist,(227)*OK+(w + 88)*OK);
Append(~Plist,(227)*OK+(w + 139)*OK);
Append(~Plist,(229)*OK+(w + 22)*OK);
Append(~Plist,(229)*OK+(w + 207)*OK);
Append(~Plist,(-3*w - 1)*OK);
Append(~Plist,(3*w - 1)*OK);
Append(~Plist,(-2*w + 19)*OK);
Append(~Plist,(2*w + 19)*OK);
Append(~Plist,(w + 17)*OK);
Append(~Plist,(w - 17)*OK);
Append(~Plist,(293)*OK+(w + 57)*OK);
Append(~Plist,(293)*OK+(w + 236)*OK);
Append(~Plist,(307)*OK+(w + 124)*OK);
Append(~Plist,(307)*OK+(w + 183)*OK);
Append(~Plist,(-6*w - 25)*OK);
Append(~Plist,(6*w - 25)*OK);
Append(~Plist,(-4*w - 27)*OK);
Append(~Plist,(4*w - 27)*OK);
Append(~Plist,(317)*OK+(w + 126)*OK);
Append(~Plist,(317)*OK+(w + 191)*OK);
Append(~Plist,(331)*OK+(w + 41)*OK);
Append(~Plist,(331)*OK+(w + 290)*OK);
Append(~Plist,(-2*w + 21)*OK);
Append(~Plist,(2*w + 21)*OK);
Append(~Plist,(349)*OK+(w + 153)*OK);
Append(~Plist,(349)*OK+(w + 196)*OK);
Append(~Plist,(-4*w + 7)*OK);
Append(~Plist,(4*w + 7)*OK);
Append(~Plist,(379)*OK+(w + 28)*OK);
Append(~Plist,(379)*OK+(w + 351)*OK);
Append(~Plist,(397)*OK+(w + 87)*OK);
Append(~Plist,(397)*OK+(w + 310)*OK);
Append(~Plist,(421)*OK+(w + 187)*OK);
Append(~Plist,(421)*OK+(w + 234)*OK);
Append(~Plist,(-7*w - 29)*OK);
Append(~Plist,(7*w - 29)*OK);
Append(~Plist,(-5*w + 33)*OK);
Append(~Plist,(-5*w - 33)*OK);
Append(~Plist,(461)*OK+(w + 165)*OK);
Append(~Plist,(461)*OK+(w + 296)*OK);
Append(~Plist,(499)*OK+(w + 32)*OK);
Append(~Plist,(499)*OK+(w + 467)*OK);
Append(~Plist,(w + 23)*OK);
Append(~Plist,(w - 23)*OK);
Append(~Plist,(509)*OK+(w + 75)*OK);
Append(~Plist,(509)*OK+(w + 434)*OK);
Append(~Plist,(-2*w + 25)*OK);
Append(~Plist,(2*w + 25)*OK);
Append(~Plist,(541)*OK+(w + 244)*OK);
Append(~Plist,(541)*OK+(w + 297)*OK);
Append(~Plist,(557)*OK+(w + 103)*OK);
Append(~Plist,(557)*OK+(w + 454)*OK);
Append(~Plist,(-5*w + 9)*OK);
Append(~Plist,(5*w + 9)*OK);
Append(~Plist,(587)*OK+(w + 253)*OK);
Append(~Plist,(587)*OK+(w + 334)*OK);
Append(~Plist,(w + 25)*OK);
Append(~Plist,(w - 25)*OK);
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

ECSearch("1.1-a",(1)*OK,[-2, 1, 1, 5, 2, 2, -6, 3, 3, 0, 0, -6, -6, 3, 3, 5, 10, 10, -12, -12, 0, 0, -16, -16, 4, 4, 15, 15, -6, -6]);
ECSearch("1.1-b",(1)*OK,[2, -1, -1, 5, -2, -2, 6, 3, 3, 0, 0, -6, -6, -3, -3, 5, -10, -10, 12, 12, 0, 0, 16, 16, 4, 4, -15, -15, -6, -6]);
ECSearch("8.1-a",(4)*OK+(2*w)*OK,[3, -1, -1, 0, 0, -2, 3, -5, 0, 8, 0, 8, 7, 3, -5, 0, 8, -8, 8, 8, 0, 0, 8, 0, -8, 11, 7, 2, 2]);
ECSearch("8.1-b",(4)*OK+(2*w)*OK,[1, -3, -1, 0, 0, 2, -5, 3, -8, 0, 8, 0, -3, -7, -5, -8, 0, -8, 8, 0, 8, -8, 0, -8, 0, -7, -11, 2, 2]);
ECSearch("8.1-c",(4)*OK+(2*w)*OK,[-1, 3, -1, 0, 0, -2, -5, 3, 8, 0, 8, 0, 3, 7, -5, 8, 0, 8, -8, 0, 8, 8, 0, -8, 0, 7, 11, 2, 2]);
ECSearch("8.1-d",(4)*OK+(2*w)*OK,[-3, 1, -1, 0, 0, 2, 3, -5, 0, -8, 0, 8, -7, -3, -5, 0, -8, 8, -8, 8, 0, 0, -8, 0, -8, -11, -7, 2, 2]);
ECSearch("9.1-a",(3)*OK,[2, 2, 2, 4, 4, -6, -6, -6, 6, 6, 0, 0, 6, 6, 14, -4, -4, 6, 6, 6, 6, -8, -8, -14, -14, -6, -6, -18, -18]);
ECSearch("9.1-b",(3)*OK,[-2, -2, -2, -4, -4, 6, -6, -6, -6, -6, 0, 0, -6, -6, 14, 4, 4, -6, -6, 6, 6, 8, 8, -14, -14, 6, 6, -18, -18]);
ECSearch("10.1-a",(w + 4)*OK,[4, 0, 3, 3, 1, 3, -2, 0, -5, 9, -1, 2, -8, 5, -10, 10, -13, -8, -10, -5, 1, -14, 4, -6, -5, 0, 14, 14]);
ECSearch("10.1-b",(w + 4)*OK,[-2, -4, 3, 1, 5, -1, -4, -4, 1, -3, -3, -4, 2, -1, 6, -4, -3, -8, 2, -11, 7, -4, 14, -20, -9, -4, 14, -6]);
ECSearch("10.1-c",(w + 4)*OK,[-4, 0, -3, -3, -1, 3, -2, 0, 5, 9, -1, -2, 8, 5, 10, -10, 13, 8, -10, -5, -1, 14, 4, -6, 5, 0, 14, 14]);
ECSearch("10.1-d",(w + 4)*OK,[2, -4, -3, -1, -5, -1, -4, 4, -1, -3, -3, 4, -2, -1, -6, 4, 3, 8, 2, -11, -7, 4, 14, -20, 9, 4, 14, -6]);
ECSearch("10.2-a",(-w + 4)*OK,[4, 0, 3, 3, 1, -2, 3, -5, 0, -1, 9, -8, 2, 5, 10, -10, -8, -13, -5, -10, -14, 1, -6, 4, 0, -5, 14, 14]);
ECSearch("10.2-b",(-w + 4)*OK,[-2, -4, 1, 3, 5, -4, -1, 1, -4, -3, -3, 2, -4, -1, -4, 6, -8, -3, -11, 2, -4, 7, -20, 14, -4, -9, -6, 14]);
ECSearch("10.2-c",(-w + 4)*OK,[-4, 0, -3, -3, -1, -2, 3, 5, 0, -1, 9, 8, -2, 5, -10, 10, 8, 13, -5, -10, 14, -1, -6, 4, 0, 5, 14, 14]);
ECSearch("10.2-d",(-w + 4)*OK,[2, -4, -1, -3, -5, -4, -1, -1, 4, -3, -3, -2, 4, -1, 4, -6, 8, 3, -11, 2, 4, -7, -20, 14, 4, 9, -6, 14]);
ECSearch("13.1-a",(13)*OK+(w)*OK,[-1, 3, -1, -5, 0, 4, -1, -1, 0, -8, -4, 8, -5, 7, -9, -4, -4, -4, -16, 0, 4, 0, -12, 16, -4, -9, 19, 18, -14]);
ECSearch("13.1-b",(13)*OK+(w)*OK,[-1, -1, 3, -5, 4, 0, -1, -1, -8, 0, 8, -4, 7, -5, -9, -4, -4, -16, -4, 4, 0, -12, 0, -4, 16, 19, -9, -14, 18]);
ECSearch("13.1-c",(13)*OK+(w)*OK,[1, 1, -3, -5, -4, 0, -1, -1, 8, 0, 8, -4, -7, 5, -9, 4, 4, 16, 4, 4, 0, 12, 0, -4, 16, -19, 9, -14, 18]);
ECSearch("13.1-d",(13)*OK+(w)*OK,[1, -3, 1, -5, 0, -4, -1, -1, 0, 8, -4, 8, 5, -7, -9, 4, 4, 4, 16, 0, 4, 0, 12, 16, -4, 9, -19, 18, -14]);
ECSearch("16.1-a",(4)*OK,[3, -1, -1, 0, 0, -2, 3, -5, 0, -8, 0, -8, 7, 3, -5, 0, -8, 8, -8, -8, 0, 0, -8, 0, 8, 11, 7, 2, 2]);
ECSearch("16.1-b",(4)*OK,[1, 1, 5, -2, -2, -6, 3, 3, 0, 0, 6, 6, 3, 3, 5, -10, -10, 12, 12, 0, 0, 16, 16, -4, -4, 15, 15, -6, -6]);
ECSearch("16.1-c",(4)*OK,[-1, 3, -1, 0, 0, -2, -5, 3, -8, 0, -8, 0, 3, 7, -5, -8, 0, -8, 8, 0, -8, -8, 0, 8, 0, 7, 11, 2, 2]);
ECSearch("16.1-d",(4)*OK,[1, -3, -1, 0, 0, 2, -5, 3, 8, 0, -8, 0, -3, -7, -5, 8, 0, 8, -8, 0, -8, 8, 0, 8, 0, -7, -11, 2, 2]);
ECSearch("16.1-e",(4)*OK,[-1, -1, 5, 2, 2, 6, 3, 3, 0, 0, 6, 6, -3, -3, 5, 10, 10, -12, -12, 0, 0, -16, -16, -4, -4, -15, -15, -6, -6]);
ECSearch("16.1-f",(4)*OK,[-3, 1, -1, 0, 0, 2, 3, -5, 0, 8, 0, -8, -7, -3, -5, 0, 8, -8, 8, -8, 0, 0, 8, 0, 8, -11, -7, 2, 2]);
ECSearch("17.1-a",(w + 3)*OK,[-1, 2, 2, -4, 4, 4, 0, 6, 6, 0, 0, -6, 12, -6, -4, -4, 8, 12, -12, 0, 12, 4, 10, 4, -14, 6, 0, 12, 6]);
ECSearch("17.1-b",(w + 3)*OK,[1, -2, -2, -4, -4, -4, 0, 6, -6, 0, 0, -6, -12, 6, -4, 4, -8, -12, 12, 0, 12, -4, -10, 4, -14, -6, 0, 12, 6]);
ECSearch("17.2-a",(w - 3)*OK,[-1, 2, 2, -4, 4, 4, 0, 6, 0, 6, -6, 0, -6, 12, -4, 8, -4, -12, 12, 12, 0, 10, 4, -14, 4, 0, 6, 6, 12]);
ECSearch("17.2-b",(w - 3)*OK,[1, -2, -2, -4, -4, -4, 0, 6, 0, -6, -6, 0, 6, -12, -4, -8, 4, 12, -12, 12, 0, -10, -4, -14, 4, 0, -6, 6, 12]);
ECSearch("20.1-a",(10)*OK+(2*w + 2)*OK,[2, 2, -4, 0, -2, 2, 2, -4, -8, 0, 0, -2, -2, 2, -8, 12, -16, 12, -8, -16, 16, -4, -8, -16, -2, 6, -6, 2]);
ECSearch("20.1-b",(10)*OK+(2*w + 2)*OK,[-2, 2, 4, 0, 2, 2, 2, 4, 8, 0, 0, 2, 2, 2, 8, -12, 16, -12, -8, -16, -16, 4, -8, -16, 2, -6, -6, 2]);
ECSearch("20.2-a",(10)*OK+(2*w + 8)*OK,[2, 2, 0, -4, -2, 2, 2, -8, -4, 0, 0, -2, -2, 2, 12, -8, 12, -16, -16, -8, -4, 16, -16, -8, 6, -2, 2, -6]);
ECSearch("20.2-b",(10)*OK+(2*w + 8)*OK,[-2, 2, 0, 4, 2, 2, 2, 8, 4, 0, 0, 2, 2, 2, -12, 8, -12, 16, -16, -8, 4, -16, -16, -8, -6, 2, 2, -6]);
ECSearch("22.1-a",(w + 2)*OK,[1, -1, 1, -2, 2, -5, -7, -4, -2, 6, -2, -1, 9, -3, -12, 4, -4, -2, -4, 2, 12, 0, 6, -12, -5, -1, 6, -18]);
ECSearch("22.1-b",(w + 2)*OK,[-1, 1, 1, 2, -2, -5, -7, 4, 2, 6, -2, 1, -9, -3, 12, -4, 4, 2, -4, 2, -12, 0, 6, -12, 5, 1, 6, -18]);
ECSearch("22.2-a",(w - 2)*OK,[-1, 1, 1, -2, 2, -7, -5, -2, -4, -2, 6, 9, -1, -3, 4, -12, -2, -4, 2, -4, 0, 12, -12, 6, -1, -5, -18, 6]);
ECSearch("22.2-b",(w - 2)*OK,[1, -1, 1, 2, -2, -7, -5, 2, 4, -2, 6, -9, 1, -3, -4, 12, 2, 4, 2, -4, 0, -12, -12, 6, 1, 5, -18, 6]);
ECSearch("23.1-a",(w + 7)*OK,[0, -3, 2, -1, -1, 2, -1, -3, 2, -4, 0, -2, 12, 6, 4, 12, 5, -14, 4, -14, 3, 0, -12, 0, -15, 4, -17, -1, 10]);
ECSearch("23.1-b",(w + 7)*OK,[2, -1, 2, -1, 1, -2, -3, -3, -6, 0, 0, 6, 0, -6, -4, -4, 11, -6, -12, -6, 15, 4, 4, 4, 13, 0, -3, -9, -6]);
ECSearch("23.1-c",(w + 7)*OK,[-2, 1, -2, -1, -1, 2, 3, -3, -6, 0, 0, 6, 0, 6, -4, 4, -11, 6, 12, -6, 15, -4, -4, 4, 13, 0, 3, -9, -6]);
ECSearch("23.1-d",(w + 7)*OK,[0, 3, -2, -1, 1, -2, 1, -3, 2, 4, 0, -2, -12, -6, 4, -12, -5, 14, -4, -14, 3, 0, 12, 0, -15, -4, 17, -1, 10]);
ECSearch("23.2-a",(w - 7)*OK,[0, 2, -3, -1, 2, -1, -1, 2, -3, 0, -4, -2, 6, 12, 4, 5, 12, 4, -14, 3, -14, -12, 0, -15, 0, -17, 4, 10, -1]);
ECSearch("23.2-b",(w - 7)*OK,[2, 2, -1, -1, -2, 1, -3, -6, -3, 0, 0, 6, -6, 0, -4, 11, -4, -12, -6, 15, -6, 4, 4, 13, 4, -3, 0, -6, -9]);
ECSearch("23.2-c",(w - 7)*OK,[-2, -2, 1, -1, 2, -1, 3, -6, -3, 0, 0, 6, 6, 0, -4, -11, 4, 12, 6, 15, -6, -4, -4, 13, 4, 3, 0, -6, -9]);
ECSearch("23.2-d",(w - 7)*OK,[0, -2, 3, -1, -2, 1, 1, 2, -3, 0, 4, -2, -6, -12, 4, -5, -12, -4, 14, 3, -14, 12, 0, -15, 0, 17, -4, 10, -1]);
ECSearch("26.1-a",(w)*OK,[1, 1, 3, 2, 2, -3, -3, -6, -6, -4, -4, -3, -3, -13, 10, 10, 2, 2, -4, -4, 0, 0, -8, -8, -19, -19, 2, 2]);
ECSearch("26.1-b",(w)*OK,[-1, -1, 3, -2, -2, -3, -3, 6, 6, -4, -4, 3, 3, -13, -10, -10, -2, -2, -4, -4, 0, 0, -8, -8, 19, 19, 2, 2]);
ECSearch("26.1-c",(w)*OK,[3, 3, -5, -6, -6, -3, -3, -2, -2, 0, 0, 7, 7, -13, 6, 6, -14, -14, 8, 8, -12, -12, -4, -4, 7, 7, -6, -6]);
ECSearch("26.1-d",(w)*OK,[-3, -3, -5, 6, 6, -3, -3, 2, 2, 0, 0, -7, -7, -13, -6, -6, 14, 14, 8, 8, 12, 12, -4, -4, -7, -7, -6, -6]);
ECSearch("32.1-a",(8)*OK+(4*w)*OK,[2, 2, -6, 0, 0, -6, 2, 2, 0, 0, 0, 0, 2, 2, -14, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -6, -6, -14, -14]);
ECSearch("32.1-b",(8)*OK+(4*w)*OK,[-2, -2, -6, 0, 0, 6, 2, 2, 0, 0, 0, 0, -2, -2, -14, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 6, 6, -14, -14]);
ECSearch("40.1-a",(2*w + 8)*OK,[4, 5, -3, 6, 2, 1, 0, -5, -6, 2, 0, -9, 2, -8, 13, -6, 4, -4, -3, -7, 4, -12, -8, 6, 2, -5, -10, -1]);
ECSearch("40.1-b",(2*w + 8)*OK,[2, 2, 6, 6, -2, -2, 6, 2, -6, -4, 0, 6, -2, -2, 2, -6, 2, -2, 12, -4, 2, 6, -8, 12, -2, -10, 14, 14]);
ECSearch("40.1-c",(2*w + 8)*OK,[-4, 5, 3, -6, -2, 1, 0, 5, 6, 2, 0, 9, -2, -8, -13, 6, -4, 4, -3, -7, -4, 12, -8, 6, -2, 5, -10, -1]);
ECSearch("40.1-d",(2*w + 8)*OK,[-2, 2, -6, -6, 2, -2, 6, -2, 6, -4, 0, -6, 2, -2, -2, 6, -2, 2, 12, -4, -2, -6, -8, 12, 2, 10, 14, 14]);
ECSearch("40.1-e",(2*w + 8)*OK,[2, -4, 3, -3, 1, 7, 0, -4, -3, -7, -3, 0, -2, -5, -10, 0, -7, -8, -6, 5, 11, 0, -14, 0, -5, -4, 14, 2]);
ECSearch("40.1-f",(2*w + 8)*OK,[-2, -4, -3, 3, -1, 7, 0, 4, 3, -7, -3, 0, 2, -5, 10, 0, 7, 8, -6, 5, -11, 0, -14, 0, 5, 4, 14, 2]);
ECSearch("40.2-a",(-2*w + 8)*OK,[4, 5, 6, -3, 2, 0, 1, -6, -5, 0, 2, 2, -9, -8, -6, 13, -4, 4, -7, -3, -12, 4, 6, -8, -5, 2, -1, -10]);
ECSearch("40.2-b",(-2*w + 8)*OK,[2, 2, 6, 6, -2, 6, -2, -6, 2, 0, -4, -2, 6, -2, -6, 2, -2, 2, -4, 12, 6, 2, 12, -8, -10, -2, 14, 14]);
ECSearch("40.2-c",(-2*w + 8)*OK,[-4, 5, -6, 3, -2, 0, 1, 6, 5, 0, 2, -2, 9, -8, 6, -13, 4, -4, -7, -3, 12, -4, 6, -8, 5, -2, -1, -10]);
ECSearch("40.2-d",(-2*w + 8)*OK,[-2, 2, -6, -6, 2, 6, -2, 6, -2, 0, -4, 2, -6, -2, 6, -2, 2, -2, -4, 12, -6, -2, 12, -8, 10, 2, 14, 14]);
ECSearch("40.2-e",(-2*w + 8)*OK,[2, -4, -3, 3, 1, 0, 7, -3, -4, -3, -7, -2, 0, -5, 0, -10, -8, -7, 5, -6, 0, 11, 0, -14, -4, -5, 2, 14]);
ECSearch("40.2-f",(-2*w + 8)*OK,[-2, -4, 3, -3, -1, 0, 7, 3, 4, -3, -7, 2, 0, -5, 0, 10, 8, 7, 5, -6, 0, -11, 0, -14, 4, 5, 2, 14]);
ECSearch("45.1-a",(15)*OK+(3*w + 3)*OK,[-2, 2, 0, 0, 2, 2, 2, 2, 2, 0, 8, 6, 6, -2, 0, 0, 10, 2, -2, 14, -12, 12, 10, -14, 18, 2, 6, 6]);
ECSearch("45.1-b",(15)*OK+(3*w + 3)*OK,[0, 2, -6, 2, 2, -6, -6, 0, 8, -8, 4, 2, 10, -6, 14, 6, -4, -8, -10, 6, 6, -14, -14, 6, 10, 10, 6, -10]);
ECSearch("45.1-c",(15)*OK+(3*w + 3)*OK,[0, -2, 6, -2, -2, -6, -6, 0, -8, -8, 4, -2, -10, -6, -14, -6, 4, 8, -10, 6, -6, 14, -14, 6, -10, -10, 6, -10]);
ECSearch("45.1-d",(15)*OK+(3*w + 3)*OK,[2, -2, 0, 0, -2, 2, 2, -2, -2, 0, 8, -6, -6, -2, 0, 0, -10, -2, -2, 14, 12, -12, 10, -14, -18, -2, 6, 6]);
ECSearch("45.2-a",(15)*OK+(3*w + 12)*OK,[-2, 2, 0, 0, 2, 2, 2, 2, 2, 8, 0, 6, 6, -2, 0, 0, 2, 10, 14, -2, 12, -12, -14, 10, 2, 18, 6, 6]);
ECSearch("45.2-b",(15)*OK+(3*w + 12)*OK,[0, 2, 2, -6, 2, -6, -6, 8, 0, 4, -8, 10, 2, -6, 6, 14, -8, -4, 6, -10, -14, 6, 6, -14, 10, 10, -10, 6]);
ECSearch("45.2-c",(15)*OK+(3*w + 12)*OK,[0, -2, -2, 6, -2, -6, -6, -8, 0, 4, -8, -10, -2, -6, -6, -14, 8, 4, 6, -10, 14, -6, 6, -14, -10, -10, -10, 6]);
ECSearch("45.2-d",(15)*OK+(3*w + 12)*OK,[2, -2, 0, 0, -2, 2, 2, -2, -2, 8, 0, -6, -6, -2, 0, 0, -2, -10, 14, -2, -12, 12, -14, 10, -2, -18, 6, 6]);
ECSearch("46.1-a",(46)*OK+(w + 16)*OK,[1, 0, 2, 1, 4, -5, 8, -3, -6, -8, 1, -9, 6, -2, 1, -12, 8, -4, 6, 4, 9, 3, 6, -18, 5, 2, 19, -4]);
ECSearch("46.1-b",(46)*OK+(w + 16)*OK,[-1, 0, 2, -1, -4, 5, 8, -3, 6, 8, 1, 9, -6, -2, -1, 12, -8, 4, 6, 4, -9, -3, 6, -18, -5, -2, 19, -4]);
ECSearch("46.2-a",(46)*OK+(w + 30)*OK,[0, 1, 2, 4, 1, -5, -3, 8, -8, -6, 1, 6, -9, -2, -12, 1, -4, 8, 4, 6, 3, 9, -18, 6, 2, 5, -4, 19]);
ECSearch("46.2-b",(46)*OK+(w + 30)*OK,[0, -1, 2, -4, -1, 5, -3, 8, 8, 6, 1, -6, 9, -2, 12, -1, 4, -8, 4, 6, -3, -9, -18, 6, -2, -5, -4, 19]);
ECSearch("49.1-a",(7)*OK,[1, -2, -2, -1, 5, -1, 0, 0, 0, -6, 6, -3, 3, 9, -3, -2, 10, -9, -3, 3, -3, 2, -10, 4, 4, 6, 6, -9, 15]);
ECSearch("49.1-b",(7)*OK,[1, -2, -2, -1, -1, 5, 0, 0, 0, 6, -6, 3, -3, -3, 9, 10, -2, -3, -9, -3, 3, -10, 2, 4, 4, 6, 6, 15, -9]);
ECSearch("49.1-c",(7)*OK,[-1, 2, 2, -1, 1, -5, 0, 0, 0, -6, 6, 3, -3, 3, -9, -10, 2, 3, 9, -3, 3, 10, -2, 4, 4, -6, -6, 15, -9]);
ECSearch("49.1-d",(7)*OK,[-1, 2, 2, -1, -5, 1, 0, 0, 0, 6, -6, -3, 3, -9, 3, 2, -10, 9, 3, 3, -3, -2, 10, 4, 4, -6, -6, -9, 15]);
ECSearch("50.2-a",(50)*OK+(w + 24)*OK,[-2, -5, 5, 0, -4, -7, 2, -1, 6, 6, 4, -7, -2, 10, -9, 4, 2, -8, 15, 5, -16, -6, -16, 6, 4, 1, 6, 9]);
ECSearch("50.2-b",(50)*OK+(w + 24)*OK,[-4, -2, 1, 3, 1, -1, -4, 4, 3, 3, 7, -2, 2, -5, -6, 8, 7, -10, -6, 11, 1, -18, 2, 0, 11, 14, -6, -12]);
ECSearch("50.2-c",(50)*OK+(w + 24)*OK,[0, -4, 3, 3, -3, 3, 6, 0, 3, -3, 3, 6, 0, 13, 6, -6, 3, 12, -10, -1, 9, 6, -4, -14, 15, -12, 6, -6]);
ECSearch("50.2-d",(50)*OK+(w + 24)*OK,[0, 1, 3, -2, 2, -7, -4, 5, -2, 2, 8, 1, -10, -12, -9, -6, 8, 12, 5, -11, 4, 16, 16, 6, 10, 13, 6, 19]);
ECSearch("50.2-e",(50)*OK+(w + 24)*OK,[-3, 5, 0, 0, -6, 3, -3, 6, -6, 6, -6, -3, -3, -5, -6, 6, -12, -12, -10, -10, 6, 6, 14, -14, -9, 9, 6, -6]);
ECSearch("50.2-f",(50)*OK+(w + 24)*OK,[-2, 0, 0, 0, 1, 3, 7, -6, 1, 6, -1, -7, -2, -10, 11, 4, -13, -8, -10, -5, -11, 4, 4, -9, 19, -9, -4, -16]);
ECSearch("50.2-g",(50)*OK+(w + 24)*OK,[0, 4, 0, 4, 5, -1, 5, 8, -5, 2, -7, 7, 2, 6, 3, 0, 5, -6, -4, -11, 7, -8, -14, -9, 19, 1, 12, 10]);
ECSearch("50.2-h",(50)*OK+(w + 24)*OK,[-2, 5, 0, -5, -4, -7, 2, -1, 1, 6, -6, 8, 8, 5, -9, -1, 12, 7, -10, 10, 4, -1, 14, 16, -6, -4, -19, 9]);
ECSearch("50.2-i",(50)*OK+(w + 24)*OK,[2, 5, 0, 5, 4, -7, 2, 1, -1, 6, -6, -8, -8, 5, 9, 1, -12, -7, -10, 10, -4, 1, 14, 16, 6, 4, -19, 9]);
ECSearch("50.2-j",(50)*OK+(w + 24)*OK,[0, 4, 0, -4, -5, -1, 5, -8, 5, 2, -7, -7, -2, 6, -3, 0, -5, 6, -4, -11, -7, 8, -14, -9, -19, -1, 12, 10]);
ECSearch("50.2-k",(50)*OK+(w + 24)*OK,[2, 0, 0, 0, -1, 3, 7, 6, -1, 6, -1, 7, 2, -10, -11, -4, 13, 8, -10, -5, 11, -4, 4, -9, -19, 9, -4, -16]);
ECSearch("50.2-l",(50)*OK+(w + 24)*OK,[3, 5, 0, 0, 6, 3, -3, -6, 6, 6, -6, 3, 3, -5, 6, -6, 12, 12, -10, -10, -6, -6, 14, -14, 9, -9, 6, -6]);
ECSearch("50.2-m",(50)*OK+(w + 24)*OK,[0, 1, -3, 2, -2, -7, -4, -5, 2, 2, 8, -1, 10, -12, 9, 6, -8, -12, 5, -11, -4, -16, 16, 6, -10, -13, 6, 19]);
ECSearch("50.2-n",(50)*OK+(w + 24)*OK,[0, -4, -3, -3, 3, 3, 6, 0, -3, -3, 3, -6, 0, 13, -6, 6, -3, -12, -10, -1, -9, -6, -4, -14, -15, 12, 6, -6]);
ECSearch("50.2-o",(50)*OK+(w + 24)*OK,[4, -2, -1, -3, -1, -1, -4, -4, -3, 3, 7, 2, -2, -5, 6, -8, -7, 10, -6, 11, -1, 18, 2, 0, -11, -14, -6, -12]);
ECSearch("50.2-p",(50)*OK+(w + 24)*OK,[2, -5, -5, 0, 4, -7, 2, 1, -6, 6, 4, 7, 2, 10, 9, -4, -2, 8, 15, 5, 16, 6, -16, 6, -4, -1, 6, 9]);
ECSearch("50.3-a",(50)*OK+(w + 26)*OK,[-2, -5, 0, 5, -4, 2, -7, 6, -1, 4, 6, -2, -7, 10, 4, -9, -8, 2, 5, 15, -6, -16, 6, -16, 1, 4, 9, 6]);
ECSearch("50.3-b",(50)*OK+(w + 26)*OK,[-4, -2, 3, 1, 1, -4, -1, 3, 4, 7, 3, 2, -2, -5, 8, -6, -10, 7, 11, -6, -18, 1, 0, 2, 14, 11, -12, -6]);
ECSearch("50.3-c",(50)*OK+(w + 26)*OK,[0, -4, 3, 3, -3, 6, 3, 3, 0, 3, -3, 0, 6, 13, -6, 6, 12, 3, -1, -10, 6, 9, -14, -4, -12, 15, -6, 6]);
ECSearch("50.3-d",(50)*OK+(w + 26)*OK,[0, 1, -2, 3, 2, -4, -7, -2, 5, 8, 2, -10, 1, -12, -6, -9, 12, 8, -11, 5, 16, 4, 6, 16, 13, 10, 19, 6]);
ECSearch("50.3-e",(50)*OK+(w + 26)*OK,[-3, 5, 0, 0, -6, -3, 3, -6, 6, -6, 6, -3, -3, -5, 6, -6, -12, -12, -10, -10, 6, 6, -14, 14, 9, -9, -6, 6]);
ECSearch("50.3-f",(50)*OK+(w + 26)*OK,[-2, 0, 0, 0, 1, 7, 3, 1, -6, -1, 6, -2, -7, -10, 4, 11, -8, -13, -5, -10, 4, -11, -9, 4, -9, 19, -16, -4]);
ECSearch("50.3-g",(50)*OK+(w + 26)*OK,[0, 4, 4, 0, 5, 5, -1, -5, 8, -7, 2, 2, 7, 6, 0, 3, -6, 5, -11, -4, -8, 7, -9, -14, 1, 19, 10, 12]);
ECSearch("50.3-h",(50)*OK+(w + 26)*OK,[-2, 5, -5, 0, -4, 2, -7, 1, -1, -6, 6, 8, 8, 5, -1, -9, 7, 12, 10, -10, -1, 4, 16, 14, -4, -6, 9, -19]);
ECSearch("50.3-i",(50)*OK+(w + 26)*OK,[2, 5, 5, 0, 4, 2, -7, -1, 1, -6, 6, -8, -8, 5, 1, 9, -7, -12, 10, -10, 1, -4, 16, 14, 4, 6, 9, -19]);
ECSearch("50.3-j",(50)*OK+(w + 26)*OK,[0, 4, -4, 0, -5, 5, -1, 5, -8, -7, 2, -2, -7, 6, 0, -3, 6, -5, -11, -4, 8, -7, -9, -14, -1, -19, 10, 12]);
ECSearch("50.3-k",(50)*OK+(w + 26)*OK,[2, 0, 0, 0, -1, 7, 3, -1, 6, -1, 6, 2, 7, -10, -4, -11, 8, 13, -5, -10, -4, 11, -9, 4, 9, -19, -16, -4]);
ECSearch("50.3-l",(50)*OK+(w + 26)*OK,[3, 5, 0, 0, 6, -3, 3, 6, -6, -6, 6, 3, 3, -5, -6, 6, 12, 12, -10, -10, -6, -6, -14, 14, -9, 9, -6, 6]);
ECSearch("50.3-m",(50)*OK+(w + 26)*OK,[0, 1, 2, -3, -2, -4, -7, 2, -5, 8, 2, 10, -1, -12, 6, 9, -12, -8, -11, 5, -16, -4, 6, 16, -13, -10, 19, 6]);
ECSearch("50.3-n",(50)*OK+(w + 26)*OK,[0, -4, -3, -3, 3, 6, 3, -3, 0, 3, -3, 0, -6, 13, 6, -6, -12, -3, -1, -10, -6, -9, -14, -4, 12, -15, -6, 6]);
ECSearch("50.3-o",(50)*OK+(w + 26)*OK,[4, -2, -3, -1, -1, -4, -1, -3, -4, 7, 3, -2, 2, -5, -8, 6, 10, -7, 11, -6, 18, -1, 0, 2, -14, -11, -12, -6]);
ECSearch("50.3-p",(50)*OK+(w + 26)*OK,[2, -5, 0, -5, 4, 2, -7, -6, 1, 4, 6, 2, 7, 10, -4, 9, 8, -2, 5, 15, 6, 16, 6, -16, -1, -4, 9, 6]);
ECSearch("52.1-a",(26)*OK+(2*w)*OK,[2, 2, -6, -2, -2, 6, 6, -6, -6, 8, 8, -6, -6, -10, -10, -10, 10, 10, -4, -4, -6, -6, -8, -8, -14, -14, 14, 14]);
ECSearch("52.1-b",(26)*OK+(2*w)*OK,[-2, -2, -6, 2, 2, 6, 6, 6, 6, 8, 8, 6, 6, -10, 10, 10, -10, -10, -4, -4, 6, 6, -8, -8, 14, 14, 14, 14]);