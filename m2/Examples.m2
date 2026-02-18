restart
load "WalshVariety.m2"
needsPackage "NumericSolutions"

P = QQ[s];

-----------------
-- Example 4.1 --
-----------------

c = (s+1)^2; factor c
BB = saturate badLocus(2,1,C=>c)
dim BB    -- 1 (projectively: 0)
degree BB -- 1
-- BB consists of the simple point corresponding to (s-1)


c = (s+1)^2*(s+2)^2; factor c
BB = badLocus(4,3,C=>c);
dim BB    -- 2 (projectively: 1)
degree BB -- 1
associatedPrimes saturate BB
-- BB consists of the projective line of degree 3 multiples of (s-1)(s-2) 
-- with two embedded points corresponding to (s-1)²(s-2) and (s-1)(s-2)²  


-----------------
-- Example 4.3 --
-----------------

c = (s+1)^3 - 1/2; d = -3*s^2 + 42*s + 21; factor(d/c) -- Simple poles -> no bad locus
I = MEP(3,2, CD=>(c,d));
dim I    -- 1 (projectively: 0)
degree I -- 7
use ring I
solveSystem(I+ideal(a_2-1))

-----------------
-- Example 5.2 --
-----------------

c = (s+2)^2*(s+3); d = (s-4)*(s-5); factor(d/c)
Mmatrix(3,1, CD => (c,d))

I = MEP(3,1, CD => (c,d))
degree I  -- 5
use ring I
solveSystem(I+ideal(a_1-1))

BB = saturate badLocus(3,1, C=>c)
dim BB    -- 1 (projectively: 0)
degree BB -- 1
-- BB consists of the simple point corresponding to (s-2)


-----------------
-- Example 5.3 --
-----------------

c = (s+2)^2*(s+3); d = (s-4)*(s-5); factor(d/c)

I = MEP(3,2, CD => (c,d));
degree I  -- 7
use ring I
solveSystem(I+ideal(a_2-1))

BB = saturate badLocus(3,2, C=>c)
dim BB    -- 1 (projectively: 0)
degree BB -- 2
-- BB consists of one double point corresponding to (s-2)