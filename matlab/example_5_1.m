% Example 5.1 (continuous time).
% Optimal first-order approximation (n = 1) of the third-order model
%     H(s) = (s^2 + 9s - 10) / (s^3 + 12s^2 + 49s + 78),
% with poles -6 and -3 +/- 2i, from Agudelo et al. (2021).
% Expected: D consists of 5 reduced points; the two real points with a
% stable denominator give the local minimizers 1.2799/(s + 9.6796), which
% is the global minimizer, and -0.0437/(s + 0.2671).
%
% Requires MacaulayLab (https://www.macaulaylab.net) on the MATLAB path.

clear; clc;

d = [1, 9, -10];
c = [1, 12, 49, 78];
sys = tf(d, c);
N = numel(c) - 1;
n = 1;

tol  = 1e-6;     % tolerance for real solutions
dEnd = 30;       % maximum degree for macaulaylab

%% Critical points
[mep, ~, M] = getMep(sys, n);

% For n = 1, M(a,c,d) is square and D is the set of roots of its
% determinant (cf. Example 3.3), here in terms of a(s) = s + a1.
disp('Roots of det M(a,c,d):')
disp(roots(sym2poly(det(M))))

sols = macaulaylab(mep, dEnd, posdim=true, verbose=true);
A = [ones(size(sols.values, 1), 1), sols.values];    % one denominator a per row
fprintf('%d critical points; Theorem 3.10 predicts %d.\n', size(A, 1), expectedDegree(N, n));

%% Reduced models for the real critical points with a stable denominator
for k = find(all(abs(imag(A)) < tol, 2)).'
    a = real(A(k, :));
    if isstable(tf(1, a))
        b = getOptimalZeros(sys, a);
        fprintf('Hhat = %.4f/(s + %.4f),  H2 error %.4f\n', b, a(2), h2Error(sys, a));
    end
end

%% H2 error over the reduced denominators a(s) = s + alpha
alpha = linspace(0.05, 12, 100);
J = arrayfun(@(x) h2Error(sys, [1, x]), alpha);

plot(alpha, J)
xlabel('a')
ylabel('H_2 error')