restart
load "WalshVariety.m2"
needsPackage "NumericSolutions"
use QQ[s]


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
-- Example 5.1 --
-----------------

d = s^2 + 9*s-10; c = s^3 + 12*s^2 + 49*s + 78; factor(d/c)
Mmatrix(3,1,CD=>(c,d))

I = MEP(3,1, CD => (c,d));
dim I     -- 1 (projectively: 0)
degree I  -- 5
solveSystem(I+ideal((ring I)_1-1))


-----------------
-- Example 5.3 --
-----------------

c = (s+2)^2*(s+3); d = (s-4)*(s-5); factor(d/c)
Mmatrix(3,1, CD => (c,d))

I = MEP(3,1, CD => (c,d));
dim I     -- 1 (projectively: 0)
degree I  -- 5
solveSystem(I+ideal((ring I)_1-1))

BB = saturate badLocus(3,1, C=>c)
dim BB                            -- 1 (projectively: 0)
degree BB                         -- 1
degree saturate(I,sub(BB,ring I)) -- 4 = 5 - 1
-- BB in DD consists of the simple point corresponding to (s-2)


-----------------
-- Example 5.4 --
-----------------

c = (s+2)^2*(s+3); d = (s-4)*(s-5); factor(d/c)

I = MEP(3,2, CD => (c,d));
dim I     -- 1 (projectively: 0)
degree I  -- 7
solveSystem(I+ideal((ring I)_2-1))

BB = saturate badLocus(3,2, C=>c)
dim BB                            -- 1 (projectively: 0)
degree BB                         -- 2
degree saturate(I,sub(BB,ring I)) -- 5 = 7 - 2
-- BB in DD consists of one double point corresponding to (s-2)


-----------------
-- Example 5.5 --
-----------------

c = (s+7)*(s+1/5); d = s + 13/11; factor(d/c)
M = Mmatrix(2,1,CD=>(c,d))

factor det M -- double root (a_0 - a_1)² and simple root (a_0 - 91/55*a_1)

Wcd = Walsh(2,1, CD=>(c,d))
decompose Wcd -- Yields the two solutions in a,b,g (both are rational)


-----------------
-- Example 5.8 --
-----------------

d = 1/256*s^3 - 9/256*s^2 + 295/256*s + 2449/256;
c = s^4 + 13*s^3 + 513/16*s^2 + 239/8*s + 153/16;
factor(d/c) -- Simple poles -> no bad locus

Walsh(4,3,CD=>(c,d))
I = MEP(4,3,CD=>(c,d));
dim I       -- 1 (projectively: 0)
degree I    -- 15
sols = solveSystem (I+ideal((ring I)_3-1))

needsPackage "RealRoots"
netList ({{"a(s)","stable"}} | for a in keys sols list (
    if not(isReal(a#0) and isReal(a#1) and isReal(a#2)) then continue;
    s' := (CC(monoid[s]))_0; as := s'^3 + a#2*s'^2 + a#1*s' + a#0;
    {as, isHurwitzStable as}
))
sub(I,matrix(QQ,{{1,3,3,1}}))
-- 3 real roots, only (1,3,3,1) is stable, and corresponds to (s+1)³
-- Note: This is not a certified computation and serves an illustrative purpose.