function J = h2Error(sys, a)
%H2ERROR  H2 error of the best reduced model with a given denominator.
%   J = h2Error(sys, a) returns ||H - b/a||_H2, where b is the optimal
%   numerator from getOptimalZeros. The stable denominator a is a row
%   vector in descending powers, as in tf.

    b = getOptimalZeros(sys, a);
    J = norm(sys - tf(b, a, sys.Ts));
end