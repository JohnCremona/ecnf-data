function NewformsAndCurves(K, bound, pbound, effort)
  primes := PrimesUpTo(pbound, K);
  allcurves := [];
  for N in IdealsUpTo(bound,K) do
        print "-------------------------------------------------";
        print "Level N of norm  ",Norm(N);
        S := HilbertCuspForms(K,N);
        newforms := NewformsOfDegree1(S);
        print #newforms, "newforms of degree 1";
        for f in newforms do
              print "-------------------------------------------------";
              print "newform: ", f;
              traces := [HeckeEigenvalue(f,P) : P in primes];
              curves := EllipticCurveSearch(N, effort : Primes:=primes, Traces:=traces);
              if #curves eq 0 then print "No curve found"; end if;
              for E in curves do
                a1,a2,a3,a4,a6:=Explode(aInvariants(E));
                printf "Curve [%o,%o,%o,%o,%o]\n",a1,a2,a3,a4,a6;
               end for;
             allcurves cat:= curves;
        end for;
  end for;
return allcurves;
end function;
