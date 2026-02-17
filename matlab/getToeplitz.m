function T = getToeplitz(a,nbRows)
    if size(a,2) ~= 1
       a = a.'; 
    end
    n = length(a)-1;
    T = toeplitz([a(1); zeros(nbRows-1,1)],[a; zeros(nbRows-1,1)]);
end