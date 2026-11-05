-- Optimality conditions from Theorem 2.2, more precisely equation (5).
-- Note that we have dehomogenized using a_n = 1 (and c_N = 1), so the equations give an irreducible variety in affine space.
Walsh = method(TypicalValue => Ideal, Options => {Discrete => false, CD => null})
Walsh (ZZ,ZZ) := o -> (N,n) -> (
    assert(N>n);
    var := (name,k) -> apply(k, i -> (getSymbol name)_i); -- variables "name" from 0 to k-1
    R := if instance(o.CD,Nothing) then (
        QQ(monoid[var("a",n),var("b",n),var("c",N),var("d",N),var("g",N-n)])
    ) else (
        (coefficientRing ring o.CD#0)(monoid[var("a",n),var("b",n),var("g",N-n)])
    );
    Rs :=  R(monoid[s]); s := Rs_0;
    L := gens R;

    a := take(L, {0, n-1});
    A := sum(n,i->a_i*s^i) + s^n;
    Atilde := if o.Discrete then sum(n,i->a_i*s^(n-i)) + 1 else A(-s);

    b := take(L, {n, 2*n-1});
    B := sum(n,i->b_i*s^i);

    (C,D) := if instance(o.CD,Nothing) then (
        c := take(L, {2*n, 2*n+N-1});
        d := take(L, {2*n+N, 2*n+2*N-1});
        (sum(N,i->c_i*s^i) + s^N, sum(N,i->d_i*s^i))
    ) else 
        (sub(o.CD#0,Rs), sub(o.CD#1,Rs));

    g := take(L, -(N-n));
    G := sum(N-n,i->g_i*s^i);

    optimalityCondition := A*D - B*C - Atilde^2*G;

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

--varvec := method(TypicalValue => Matrix)
--varvec (IndexedVariableTable) := y -> (
--    assert instance(y, IndexedVariableTable);
--    vs := select(values y, v -> instance(class v,PolynomialRing));
--    return matrix apply(vs, v -> {v});
--)

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
    var := (name,k) -> apply(k, i -> (getSymbol name)_i); -- variables "name" from 0 to k-1
    R := if instance(o.CD,Nothing) then (
        QQ(monoid[var("a",n+1),var("c",N+1),var("d",N)])
    ) else (
        (coefficientRing ring o.CD#0)(monoid[var("a",n+1)])
    );

    aa := transpose (vars R)_{0..n};
    atilde := if o.Discrete then
        matrix(for i to n list {aa_(n-i,0)})
    else
        diagonalMatrix(for i to n list (-1)^i) * aa;

    (cc,dd) := if instance(o.CD,Nothing) then
        (transpose (vars R)_{n+1..n+N+1}, transpose (vars R)_{n+N+2..n+2*N+1})
    else
        (coeffvec(o.CD#0), coeffvec(o.CD#1));

    return matrix entries ( (transpose T(aa,N) * dd) | (transpose T(cc,n)) | (transpose T(atilde,N) * transpose T(atilde,N-n)) );
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
    var := (name,k) -> apply(k, i -> (getSymbol name)_i); -- variables "name" from 0 to k-1
    R := if instance(o.C,Nothing) then (
        QQ(monoid[var("a",n+1),var("c",N+1)])
    ) else (
        (coefficientRing ring o.C)(monoid[var("a",n+1)])
    );

    aa := transpose (vars R)_{0..n};
    aatilde := if o.Discrete then
        matrix(for i to n list {a_(n-i)})
    else
        diagonalMatrix(for i to n list (-1)^i) * aa;
    cc := if instance(o.C,Nothing) then
        transpose (vars R)_{n+1..N+1}
    else
        coeffvec(o.C);

    Bmatrix := (transpose T(cc,n)) | (transpose T(aatilde,N) * transpose T(aatilde,N-n));
    return minors(N, Bmatrix);
)

-- The count of critical points from Theorem 3.10
expectedDegree = method(TypicalValue => ZZ)
expectedDegree (ZZ,ZZ) := (N,n) -> sum(for j from 0 to n list binomial(N-n-1+j,j)*2^j)
