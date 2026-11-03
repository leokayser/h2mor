% Example 5.5 (critical point with multiplicity).
% First-order approximation (n = 1) of
%     H(s) = (11s + 13) / ((s + 7)(s + 1/5)).
% Expected: det M(a,c,d) has a double root at a = s + 1 and a simple root
% at a = s + 91/55, so D has 2 points while deg D = 3. Both critical points
% are real and stable.
%
% Requires MacaulayLab (https://www.macaulaylab.net) on the MATLAB path.

clear; clc;

d = [11, 13];
c = poly([-7, -1/5]);
sys = tf(d, c);
N = numel(c) - 1;
n = 1;

dEnd = 30;       % maximum degree for macaulaylab

%% Critical points
[mep, ~, M] = getMep(sys, n);

% For n = 1, D is the set of roots of det M(a,c,d), here in terms of
% a(s) = s + a1; the root a1 = 1 appears twice.
disp('Roots of det M(a,c,d):')
disp(roots(sym2poly(det(M))))

sols = macaulaylab(mep, dEnd, posdim=true, verbose=true);
disp('Points of D, as a1 in a(s) = s + a1:')
disp(sols.values)
fprintf('deg D = %d (Theorem 3.10).\n', expectedDegree(N, n));

%% H2 error over the reduced denominators a(s) = s + alpha
alpha = linspace(0.5, 2, 50);
J = arrayfun(@(x) h2Error(sys, [1, x]), alpha);

plot(alpha, J)
xlabel('\alpha')
ylabel('H_2 error')