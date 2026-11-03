% Example 5.3 (non-empty bad locus).
% First-order approximation (n = 1) of
%     H(s) = (s - 4)(s - 5) / ((s + 2)^2 (s + 3)),
% whose denominator has a double root.
% Expected: D consists of 5 points. One of them, a = s - 2, forms the bad
% locus B, which leaves 4 critical points (Proposition 4.4).
%
% Requires MacaulayLab (https://www.macaulaylab.net) on the MATLAB path.

clear; clc;

d = poly([4, 5]);
c = poly([-2, -2, -3]);
sys = tf(d, c);
N = numel(c) - 1;
n = 1;

dEnd = 30;       % maximum degree for macaulaylab

%% Determinantal variety D
[mep, ~, M] = getMep(sys, n);
sols = macaulaylab(mep, dEnd, posdim=true, verbose=true);
disp('Points of D, as a1 in a(s) = s + a1:')
disp(sols.values)

%% Bad locus B: points where B(a,c), the last N columns of M(a,c,d), drops rank
% A null vector of B(a,c) is normalized as (1; q).
B = M(:, 2:end);
q = sym('q', [N-1, 1]);
mepBad = symToMep(B * [1; q], q);
solsBad = macaulaylab(mepBad, dEnd, posdim=true, verbose=true);
disp('Points of B:')
disp(solsBad.values)

%% Counts
fprintf('Generic count (14): %d; upper bound from Proposition 4.4: %d.\n', ...
    expectedDegree(N, n), expectedDegree(N, n) - nchoosek(N-2, n-1));