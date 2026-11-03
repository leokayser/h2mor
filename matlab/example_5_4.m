% Example 5.4 (multiplicity in the bad locus points).
% Second-order approximation (n = 2) of the model of Example 5.3,
%     H(s) = (s - 4)(s - 5) / ((s + 2)^2 (s + 3)).
% Expected: D consists of 6 points, while deg D = 7. The point
% a = s^2 - 5s + 6 has multiplicity two and forms the bad locus B, which
% leaves 5 critical points (Corollary 4.3); Proposition 4.4 only gives the
% upper bound 6.
%
% Requires MacaulayLab (https://www.macaulaylab.net) on the MATLAB path.

clear; clc;

d = poly([4, 5]);
c = poly([-2, -2, -3]);
sys = tf(d, c);
N = numel(c) - 1;
n = 2;

dEnd = 30;       % maximum degree for macaulaylab

%% Determinantal variety D
[mep, ~, M] = getMep(sys, n);
sols = macaulaylab(mep, dEnd, posdim=true, verbose=true);
disp('Points of D, as (a1, a2) in a(s) = s^2 + a1 s + a2:')
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