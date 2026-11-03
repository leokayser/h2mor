function m = expectedDegree(N, n)
%EXPECTEDDEGREE  Number of critical points for generic data (Theorem 3.10).
%   m = expectedDegree(N, n) evaluates (14) for a higher-order model of
%   McMillan degree N and reduced order n.

    j = 0:n;
    m = sum(arrayfun(@(k) nchoosek(N-n-1+k, k), j) .* 2.^j);
end