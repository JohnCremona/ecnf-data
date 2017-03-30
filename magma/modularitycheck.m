

// This is based on theorems in the following papers
// 
// [Chen] Imin Chen, The Jacobian of Modular Curves Associated
// to Cartan Subgroups, Oxford DPhil thesis, 1996.
//
// [FLS] Nuno Freitas, Bao Le Hung, and S. Siksek, 
// Elliptic curves over Real Quadratic Fields are Modular,
// Invent. math. (2015) 201, pp. 159--206. 
//
// [Thorne] Jack Thorne, Automorphy of some residually dihedral Galois representations. 
// Mathematische Annalen 364 (2016), No. 1--2, pp. 589--648 
//


// Given an elliptic curve $E$ over a totally real field $K$
// this returns either true or false. 
// true means that E is proven to be modular.
// false means that E might not yet be proven to be modular.

isModular:=function(E);
	K:=BaseRing(E);
	assert IsTotallyReal(K);
	if Degree(K) in {1,2} then
		return true;
	end if;
	if HasComplexMultiplication(E) then
		return true;
	end if;
	j:=jInvariant(E);
	Kx<x>:=PolynomialRing(K);
	//
	// We look first at the mod 3 represenation.
	//
	c3:=Evaluate(ClassicalModularPolynomial(3),[x,j]);
	//
	// The mod 3 representation of E is reducible if and only if c3 has a root.
	//
	s3:=(x-9)^3*(x+3)^3-j*x^3;
	// 
	// The image of the mod 3 representation lies in the normalizer of a split Cartan
	// subgroup if and only if s3 has a root (see Chen's thesis, page 68).
	//
	if #Roots(c3) eq 0 and #Roots(s3) eq 0	
		// Modularity now follows from [FLS] Theorem 3 and Proposition 4.1.
		then return true;
	end if;
	// 
	// Next we look at the mod 5 representation.
	//
	c5:=Evaluate(ClassicalModularPolynomial(5),[x,j]);
	//
	// The mod 5 represenation is irreducible if and only if c5 has a root.
	//
	//
	ns5:=5^4*(2*x+1)*(x+1)^3*(6*x^2+21*x+19)^3-j*(x^2+x-1)^5;
	//
	// The image of the mod 5 represenation is contained
	// in the normalizer of a non-split Cartan subgroup
	// if and only if ns5 has a root [Chen, p68].
	//
	s5:=((x^2-5)*(x^2+5*x+10)*(x+5))^3-j*(x^2+5*x+5)^5;
	//
	// The image of the mod 5 represenation is contained
	// in the normalizer of a split Cartan subgroup
	// if and only if s5 has a root [Chen, p68].
	//
	if #Roots(c5) eq 0 then
		if #Roots(x^2-5) eq 0 then
			// \sqrt{5} does not belong to K
			// Modularity now follows from [Thorne] Theorem 1.1
			return true;
		else
			//
			// \sqrt{5} belongs to K and Thorne's Theorem is 
			// inapplicable. 
			if #Roots(ns5) eq 0 and #Roots(s5) eq 0 then
				// Modularity follows from [FLS] Proposition 2.1 
				// and Theorem 3. 
				return true;
			end if;
		end if;	
	end if;
	// Ditto with the mod 7 representation.		
	c7:=Evaluate(ClassicalModularPolynomial(7),[x,j]);
	ns7:=((4*x^2+5*x+2)*(x^2+3*x+4)*(x^2+10*x+4)*(3*x+1))^3-j*(x^3+x^2-2*x-1)^7;
	s7:=((x^2-5*x+8)*(x^2-5*x+1)*(x^4-5*x^3+8*x^2-7*x+7)*(x+1))^3*x-j*(x^3-4*x^2+3*x+1)^7;
	//
	if #Roots(c7) eq 0 and #Roots(ns7) eq 0 and #Roots(s7) eq 0 then
		return true;
		// FLS Theorem 4 and Lemma 2.2.
	else
		return false; // We've run out of tricks!
	end if;
end function;


/*

// Uncomment to run a sanity check.
// The purpose of this check is convince ourselves that the formulae
// for s3, ns5, s5, ns7, s7 are sensible (and have been copied correctly
// from Chen's thesis!). We make use of the fact that an elliptic curve
// with complex multiplication by an order R has mod p image contained
// in the non-split Cartan normalizer if p is inert in R and 
// in the split Cartan normalizer if p is split in R. 

//
//

CM:=[<-3, 3, -12288000>, <-3, 2, 54000>, <-3, 1, 0>, <-4, 2, 287496>, 
	<-4, 1, 1728>, <-7, 2, 16581375>, <-7, 1, -3375>, 
	<-8, 1, 8000>, <-11, 1, -32768>, <-19, 1, -884736>, 
	<-43, 1, -884736000>, <-67, 1, -147197952000>, 
	<-163, 1, -262537412640768000>
];
// This a list <D,f,j> where D is fundamental CM discriminant, f is conductor (what's that?)
// and j is jInvariant. This gives the 13 CM j-invariants over \Q

Qx<x>:=PolynomialRing(Rationals());
for trip in CM do
	D,_,j:=Explode(trip);
	if KroneckerSymbol(D,3) eq 1 then
		s3:=(x-9)^3*(x+3)^3-j*x^3;	
		assert #Roots(s3) ge 1;
	end if;
	if KroneckerSymbol(D,5) eq 1 then
		s5:=((x^2-5)*(x^2+5*x+10)*(x+5))^3-j*(x^2+5*x+5)^5;
		assert #Roots(s5) ge 1;
	end if;
	if KroneckerSymbol(D,7) eq 1 then
		s7:=((x^2-5*x+8)*(x^2-5*x+1)*(x^4-5*x^3+8*x^2-7*x+7)*(x+1))^3*x-j*(x^3-4*x^2+3*x+1)^7;
		assert #Roots(s7) ge 1;
	end if;
	if KroneckerSymbol(D,5) eq -1 then
		ns5:=5^4*(2*x+1)*(x+1)^3*(6*x^2+21*x+19)^3-j*(x^2+x-1)^5;
		assert #Roots(ns5) ge 1;
	end if;
	if KroneckerSymbol(D,7) eq -1 and j ne 1728 then
		ns7:=((4*x^2+5*x+2)*(x^2+3*x+4)*(x^2+10*x+4)*(3*x+1))^3-j*(x^3+x^2-2*x-1)^7;
		assert #Roots(ns7) ge 1;
	end if;
end for;

// Reason for excluding j=1728 in the last check:
// ns7 has the form
// (1728-j)*x^21+l.o.t
// so the modular curve ns7=0 has a rational point with j=1728 and x=infinity.

*/
