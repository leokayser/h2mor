-- Optimality conditions from Theorem 2.2.
-- Note that we have dehomogenized using a_n = 1 (and c_N = 1), so the equations give an irreducible variety in affine space.
Walsh = method(TypicalValue => Ideal, Options => {Discrete => false})
Walsh (ZZ,ZZ) := o -> (N,n) -> (
    assert(N>n);
    R := QQ[a_0..a_(n-1),b_0..b_(n-1),c_0..c_(N-1),d_0..d_(N-1),g_0..g_(N-n-1)];
    Rs :=  R[s];

    A := sum(n,i->a_i*s^i) + s^n;
    B := sum(n,i->b_i*s^i);
    C := sum(N,i->c_i*s^i) + s^N;
    D := sum(N,i->d_i*s^i);
    G := sum(N-n,i->g_i*s^i);

    optimalityCondition = if o.Discrete then
        A*D-B*C-(A(-s))^2*G
    else
        A*D-B*C-(numerator(A(1/s)))^2*G;

    return sub(ideal last coefficients optimalityCondition, R);
)

vecToRingList := method(TypicalValue => Sequence)
vecToRingList (Matrix) := yy -> (ring yy, flatten entries yy)

T := method(TypicalValue => Matrix) 
T (Matrix,ZZ) := (aa,k) -> ( -- Expects aa to be a vector (= 1x(r+1) or (r+1)x1 matrix)
    assert (k >= 0);
    (R,a) := vecToRingList(aa);
    r := #a - 1;
    Toep := mutableMatrix(R, k, r+k);
    for i from 0 to k-1 do for j from 0 to r do Toep_(i,i+j) = a#j;
    return matrix Toep;
)

varvec := method(TypicalValue => Matrix)
varvec (IndexedVariableTable) := y -> (
    assert instance(y, IndexedVariableTable);
    vs := select(values y, v -> instance(class v,PolynomialRing));
    return matrix apply(vs, v -> {v});
)

coeffvec := method(TypicalValue => Matrix)
coeffvec (RingElement) := f -> (
    t := (ring f)_0;
    return matrix(for i to (degree f)#0 list {coefficient(t^i,f)});
)

-- Matrix M(a;c,d) from section 3.3
-- Option: Discrete (false (default) or true): Switches between discrete or continuous time
--         CD (null (default) or pair (c,d) of polynomials): Whether to use a specific data point H = d/c
Mmatrix = method(TypicalValue => Ideal, Options => {Discrete => false, CD => null})
Mmatrix (ZZ,ZZ) := o -> (N,n) -> (
    assert(N>n);

    R := if instance(o.CD,Nothing) then
        QQ(monoid[a_0..a_n,c_0..c_N,d_0..d_(N-1)])
    else
        QQ(monoid[a_0..a_n]);

    aa := transpose (vars R)_{0..n};
    atilde := if o.Discrete then
        matrix(for i to n list {aa_(n-i,0)})
    else
        diagonalMatrix(for i to n list (-1)^i) * aa;

    (cc,dd) := if instance(o.CD,Nothing) then
        (varvec c, varvec d)
    else
        (coeffvec(o.CD#0), coeffvec(o.CD#1));

    return (transpose T(aa,N) * dd) | (transpose T(cc,n)) | (transpose T(atilde,N) * transpose T(atilde,N-n));
)

-- The MEP from section 3.3 given as the ideal of minors of M(a;c,d)
-- Option: Discrete (false (default) or true): Switches between discrete or continuous time
--         CD (null (default) or pair (c,d) of polynomials): Whether to use a specific data point H = d/c
MEP = method(TypicalValue => Ideal, Options => {Discrete => false, CD => null})
MEP (ZZ,ZZ) := o -> (N,n) -> minors(N+1, Mmatrix(N,n, Discrete=>o.Discrete, CD=>o.CD)) ;

-- The bad locus from section 3.3 given as the ideal of minor of B(a;c)
-- Option: Discrete (false (default) or true): Switches between discrete or continuous time
--         C (null (default) or polynomial c): Whether to use a specific data point c
badLocus = method(TypicalValue => Ideal, Options => {Discrete => false, C => null})
badLocus (ZZ,ZZ) := o -> (N,n) -> (
    assert(N>n);

    R := if instance(o.C,Nothing) then
        QQ(monoid[a_0..a_n,c_0..c_N])
    else
        QQ(monoid[a_0..a_n]);

    aa := transpose (vars R)_{0..n};
    aatilde := if o.Discrete then
        matrix(for i to n list {a_(n-i)})
    else
        diagonalMatrix(for i to n list (-1)^i) * aa;
    cc := if instance(o.C,Nothing) then
        varvec c
    else
        coeffvec(o.C);

    Bmatrix := (transpose T(cc,n)) | (transpose T(aatilde,N) * transpose T(aatilde,N-n));
    return minors(N, Bmatrix);
)

-- The count of critical points from Theorem 3.6
expectedDegree = method(TypicalValue => ZZ)
expectedDegree (ZZ,ZZ) := (N,n) -> sum(for j from 0 to n list binomial(N-n-1+j,j)*2^j)
