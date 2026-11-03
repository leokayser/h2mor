function b = getOptimalZeros(sys, a)
%GETOPTIMALZEROS  Optimal numerator for a fixed reduced denominator.
%   b = getOptimalZeros(sys, a) returns the numerator b that minimizes the
%   H2 error ||H - b/a|| for the given stable denominator a, by solving the
%   linear system (21),
%       [T_n(c)'  T_N(atilde)'] [b; q] = T_N(a)' d.
%   Its coefficient matrix is the Sylvester matrix of c and atilde. The
%   polynomials a and b are row vectors in descending powers, as in tf.

    [c, d] = tfCoefficients(sys);
    N = numel(d);
    a = flip(a(:));                       % ascending order, as in the paper
    n = numel(a) - 1;

    if sys.Ts == 0
        atilde = a .* (-1).^(0:n).';      % atilde(s) = a(-s)
    else
        atilde = flip(a);                 % atilde(z) = z^n a(1/z)
    end

    x = [getToeplitz(c, n).', getToeplitz(atilde, N).'] \ (getToeplitz(a, N).' * d);
    b = flip(x(1:n)).';
end