% Example 5.6 (prescribed multiplicity structure).
% Model (20) of Example 5.2, with the reduced model restricted to a double
% pole, a(z) = (z - p)^2, i.e., to the discriminant curve in Figure 1.
% On this curve the stationarity conditions become (22), a polynomial
% eigenvalue problem M'(p) of size 8.
% Expected: det M'(p) has 17 roots, of which 7 are real and 3 lie in the
% unit disk, p = -0.6196, 0.1478, 0.7503 (the local minimizers on the curve).
% The local maximizers on the curve, p = -0.4556 and 0.5974, are real stable
% critical points of the unrestricted problem with n = 1.
%
% Requires MacaulayLab (https://www.macaulaylab.net) on the MATLAB path.

clear; clc;

d = [0.0448, 0.2368, 0.0013, 0.0211, 0.2250, 0.0219];
c = [1, -1.2024, 2.3675, -2.0039, 2.2337, -1.0420, 0.8513];
sys = tf(d, c, -1);

tol = 1e-6;      % tolerance for real solutions

%% Restricted problem: a(z) = (z - p)^2
mep = getRestrictedMep(sys, 2);
sols = macaulaylab(mep, 50, posdim=true, verbose=true);

p = sols.values;
pReal = real(p(abs(imag(p)) < tol));
pStable = pReal(abs(pReal) < 1);
fprintf('Restricted problem: %d solutions, %d real, %d in the unit disk.\n', ...
    numel(p), numel(pReal), numel(pStable));

% Coordinates (a1, a0) in Figure 1, since (z - p)^2 = z^2 - 2p z + p^2
disp('Local minimizers on the curve, as (p, a1, a0):')
disp([pStable, -2*pStable, pStable.^2])

%% Unrestricted problem with n = 1
% The reduced model b/(z + a1) has its pole at -a1.
mep1 = getMep(sys, 1);
sols1 = macaulaylab(mep1, 30, posdim=true, verbose=true);

poles = -sols1.values;
poles = real(poles(abs(imag(poles)) < tol));
disp('Real stable critical points for n = 1 (local maximizers on the curve):')
disp(poles(abs(poles) < 1))