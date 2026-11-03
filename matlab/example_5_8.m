% Example 5.8 (optimal reduced model with a triple pole).
% Third-order approximation (n = 3) of the fourth-order model of Van Dooren
% et al. (2010, Section 6), given there in state-space form,
%     H(s) = (s^3 - 9s^2 + 295s + 2449) / (256 (s^4 + 13s^3 + 513/16 s^2 + 239/8 s + 153/16)).
% Expected: D consists of 15 points. Three of them are real, and only
% a = (s + 1)^3 is stable; it gives the global minimizer Hhat = 1/(s + 1)^3.
% Restricted to a(s) = (s - p)^3, det M'(p) has 15 roots, of which 5 are
% real and 3 stable, p = -1, -0.084, -0.476.
%
% Requires MacaulayLab (https://www.macaulaylab.net) on the MATLAB path.

clear; clc;

d = [1, -9, 295, 2449] / 256;
c = [1, 13, 513/16, 239/8, 153/16];
sys = tf(d, c);
N = numel(c) - 1;
n = 3;

tol  = 1e-6;     % tolerance for real solutions
dEnd = 30;       % maximum degree for macaulaylab

%% Critical points
mep = getMep(sys, n);
sols = macaulaylab(mep, dEnd, posdim=true, verbose=true);
A = [ones(size(sols.values, 1), 1), sols.values];    % one denominator a per row

isReal = all(abs(imag(A)) < tol, 2);
fprintf('%d critical points (Theorem 3.10: %d), %d real.\n', ...
    size(A, 1), expectedDegree(N, n), nnz(isReal));

for k = find(isReal).'
    a = real(A(k, :));
    if isstable(tf(1, a))
        b = getOptimalZeros(sys, a);
        fprintf('Stable: a = [%s], b = [%s], H2 error %.4f\n', ...
            num2str(a, ' %.4f'), num2str(b, ' %.4f'), h2Error(sys, a));
    end
end

%% Restricted problem: a(s) = (s - p)^3
mepR = getRestrictedMep(sys, n);
solsR = macaulaylab(mepR, 50, posdim=true, verbose=true);

p = solsR.values;
pReal = real(p(abs(imag(p)) < tol));
fprintf('Restricted problem: %d solutions, %d real.\n', numel(p), numel(pReal));
disp('Stable solutions p:')
disp(pReal(pReal < 0))