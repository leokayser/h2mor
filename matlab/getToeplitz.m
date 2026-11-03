function T = getToeplitz(f, l)
%GETTOEPLITZ  Banded Toeplitz matrix T_l(f) from the notation section of the paper.
%   T = getToeplitz(f, l) returns the l-by-(l+k) matrix whose rows are l
%   successive shifts of the coefficient vector f = (f0, ..., fk), given in
%   ascending order. For a polynomial h, T_{k+1}(h).' * f is the coefficient
%   vector of the product h*f. Works for numeric and symbolic f.

    f = f(:).';
    k = numel(f) - 1;
    T = repmat(0*f(1), l, l + k);
    for i = 1:l
        T(i, i:i+k) = f;
    end
end