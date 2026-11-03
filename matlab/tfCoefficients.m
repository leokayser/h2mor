function [c, d] = tfCoefficients(sys)
%TFCOEFFICIENTS  Coefficient vectors of H = d/c in ascending order, with c monic.
%   [c, d] = tfCoefficients(sys) returns, for a strictly proper SISO tf
%   object sys of order N, the column vectors c (N+1 entries) and d (N
%   entries) such that c(s) = c0 + c1 s + ... + s^N and
%   d(s) = d0 + d1 s + ... + d_{N-1} s^(N-1).

    [num, den] = tfdata(sys, 'v');
    N = numel(den) - 1;
    num = [zeros(1, N + 1 - numel(num)), num];     % pad to length N+1
    c = flip(den).' / den(1);
    d = flip(num(2:end)).' / den(1);
end