function b = getOptimalZeros(c,d,a,isContinuous)
    n = length(a)-1;
    N = length(c)-1;

    Td = getToeplitz(d,n+1);
    Tc = getToeplitz(c,n);

    if isContinuous
        a_t = zeros(size(a));
        for i = 1:length(a)
            a_t(i) = a(i)*(-1)^(i-1);
        end
        Ta_t = getToeplitz(a_t,N);
    else
        Ta_t = getToeplitz(flip(a),N);
    end

    tmp = [Tc.', Ta_t.']\(Td.'*a);
    b = tmp(1:n);
end
