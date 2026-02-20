restart
load "WalshVariety.m2"
needsPackage "NumericalCertification"

use QQ[s]
c = s^4 +  13 *s^3 + 513/16 *s^2 + 239/8  *s +  153/16;
d =      1/256*s^3 -   9/256*s^2 + 295/256*s + 2449/256;
factor(d/c) -- Simple poles -> no bad locus

I = MEP(4,3, CD=>(c,d));
dim I       -- 1 (projectively: 0)
degree I    -- 15
sols = solveSystem(I+ideal((ring I)_3-1))

needsPackage "RealRoots"
for a in keys sols do (
    if not(isReal(a#0) and isReal(a#1) and isReal(a#2)) then continue;
    s' := (CC(monoid[s]))_0; as := s'^3 + a#2*s'^2 + a#1*s' + a#0;
    print(as, isHurwitzStable as);
)
sub(I,matrix(QQ,{{1,3,3,1}}))
-- 3 real roots, only (1,3,3,1) is stable, and corresponds to (s+1)³
-- Note: This is not a certified computation and serves an illustrative purposes.
-- A certified computation can be obtained