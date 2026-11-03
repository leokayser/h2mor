% Example 5.2 (discrete time).
% Second-order approximation (n = 2) of the sixth-order discrete-time model
% (20), from Lagauw et al. (2023, Example 2).
% Expected: 49 critical points, as predicted by Theorem 3.10, of which 11 are
% real and 5 correspond to stable models. The second part computes the
% contour plot of Figure 1.
%
% Requires MacaulayLab (https://www.macaulaylab.net) on the MATLAB path.

clear; clc;

d = [0.0448, 0.2368, 0.0013, 0.0211, 0.2250, 0.0219];
c = [1, -1.2024, 2.3675, -2.0039, 2.2337, -1.0420, 0.8513];
sys = tf(d, c, -1);
N = numel(c) - 1;
n = 2;

tol  = 1e-6;     % tolerance for real solutions
dEnd = 30;       % maximum degree for macaulaylab

%% Critical points
mep = getMep(sys, n);
sols = macaulaylab(mep, dEnd, posdim=true, verbose=true);
A = [ones(size(sols.values, 1), 1), sols.values];    % one denominator a per row

isReal = all(abs(imag(A)) < tol, 2);
isStab = false(size(isReal));
for k = find(isReal).'
    isStab(k) = isstable(tf(1, real(A(k, :)), -1));
end
fprintf('%d critical points (Theorem 3.10: %d), %d real, %d stable.\n', ...
    size(A, 1), expectedDegree(N, n), nnz(isReal), nnz(isStab));

%% Contour plot of Figure 1
% H2 error over the stable reduced denominators a(z) = z^2 + a1 z + a0,
% each with its optimal numerator b.
nGrid = 50;
margin = 1e-6;
a1 = linspace(-2, 2, nGrid);
a0 = linspace(-1, 1, nGrid);

J = nan(nGrid);
for i = 1:nGrid
    for j = 1:nGrid
        a = [1, a1(j), a0(i)];
        if max(abs(roots(a))) < 1 - margin
            J(i, j) = h2Error(sys, a);
        end
    end
end
fprintf('H2 norm of H: %.4f\n', norm(sys));

levels = [0.86, 0.88, 0.9, 0.94, 0.98, 1, 1.04, 1.08, 1.1, 1.12, 1.16, ...
          1.2, 1.24, 1.26, 1.27, 1.28, 1.285, 1.289, 1.291, 1.2915, 1.2918];
contour(a1, a0, J, levels)
colormap jet
colorbar
xlabel('a_1')
ylabel('a_0')