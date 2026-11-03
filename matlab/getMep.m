function [mep, eqs, M] = getMep(sys, n)
%GETMEP  MEP formulation of the critical equations of H2 model reduction.
%   [mep, eqs, M] = getMep(sys, n) returns, for a SISO tf object sys and a
%   reduced order n, the critical equations (6)
%       T_N(a)' d - T_n(c)' b - T_N(atilde)' T_{N-n}(atilde)' g = 0
%   as the symbolic vector eqs, the matrix M(a,c,d) of (9), and the
%   corresponding MacaulayLab MEP. The model is continuous-time if
%   sys.Ts == 0 and discrete-time otherwise (Section 4.2).
%
%   The reduced denominator is monic, a(s) = s^n + a1 s^(n-1) + ... + an,
%   so each eigentuple (a1, ..., an) gives the denominator [1, a1, ..., an]
%   in the descending order used by tf. The bad locus is described by
%   B(a,c) = M(:, 2:end).

    [c, d] = tfCoefficients(sys);
    N = numel(d);
    isContinuous = (sys.Ts == 0);

    % Unknowns. a and b are indexed in descending order, as in tf; the
    % vectors below hold their coefficients in ascending order, as in the paper.
    a = flip([1, sym('a', [1 n])]).';
    b = flip(sym('b', [1 n])).';
    g = sym('g', [N-n, 1]);

    if isContinuous
        atilde = a .* (-1).^(0:n).';    % atilde(s) = a(-s)
    else
        atilde = flip(a);               % atilde(z) = z^n a(1/z)
    end

    % Critical equations (6) and the matrix M(a,c,d) of (9)
    M = [getToeplitz(a, N).' * d, ...
         getToeplitz(c, n).', ...
         getToeplitz(atilde, N).' * getToeplitz(atilde, N-n).'];
    eqs = M * [1; -b; -g];

    mep = symToMep(eqs, [b; g]);
    M = vpa(M);
end