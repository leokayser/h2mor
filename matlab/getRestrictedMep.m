function [mep, eqs, M] = getRestrictedMep(sys, n)
%GETRESTRICTEDMEP  Critical equations restricted to reduced models with a single pole.
%   [mep, eqs, M] = getRestrictedMep(sys, n) restricts the reduced
%   denominator to a = (s - p)^n, i.e., a single pole p of multiplicity n,
%   and replaces atilde^2 in (5) by
%       atilde' = (s + p)^(n+1)      in continuous time,
%       atilde' = (1 - p z)^(n+1)    in discrete time,
%   as in Examples 5.6 and 5.8. The critical equations become
%       T_N(a)' d - T_n(c)' b - T_{N-1}(atilde')' h = 0,
%   and M = M'(p) is square of size N + n, so that the problem is a
%   polynomial eigenvalue problem in p. As in getMep, b is indexed in
%   descending order.

    [c, d] = tfCoefficients(sys);
    N = numel(d);
    isContinuous = (sys.Ts == 0);

    p = sym('p');
    b = flip(sym('b', [1 n])).';
    h = sym('h', [N-1, 1]);

    % Coefficients in ascending order
    a = (binomials(n) .* (-p).^(n:-1:0)).';                  % (s - p)^n
    if isContinuous
        atilde = (binomials(n+1) .* p.^(n+1:-1:0)).';        % (s + p)^(n+1)
    else
        atilde = (binomials(n+1) .* (-p).^(0:n+1)).';        % (1 - p z)^(n+1)
    end

    M = [getToeplitz(a, N).' * d, ...
         getToeplitz(c, n).', ...
         getToeplitz(atilde, N-1).'];
    eqs = M * [1; -b; -h];

    mep = symToMep(eqs, [b; h]);
    M = vpa(M);
end

function v = binomials(k)
% Binomial coefficients nchoosek(k, 0), ..., nchoosek(k, k)
    v = arrayfun(@(j) nchoosek(k, j), 0:k);
end