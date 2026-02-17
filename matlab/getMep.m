function [mep,p,M] = getMep(sys,m)
    
    % Parse input:
    c = sys.Denominator{1}(2:end);
    d = sys.Numerator{1}(2:end);
    n = length(d);
    isContinuous = (sys.Ts == 0);

    % Compose Ta coefficient matrix:
    a_sym = sym('a',[1, m]); % Symbolic coefficients
    a = [1,a_sym];
    Ta = getToeplitz(flip(a), n);

   
    % Compose Tc coefficient matrix:
    tmp = [1, c];
    Tc = getToeplitz(flip(tmp), m);

    % Compose Ta_tilde matrix:
    if isContinuous
        tmp = flip(a);
        for i = 1:m+1
            tmp(i) = (-1)^(i-1) * tmp(i);
        end
    else
        tmp = a;
    end
    Tat = getToeplitz(tmp, n);
    Tatt = getToeplitz(tmp, n-m);
    
    syms b [1 m]
    syms g [1 n-m]
    
    % Get system of equations:
    p = (Ta).'*flip(d).' - (Tc).'*flip(b).' - (Tat).'*(Tatt).'*g.';
    M = vpa([(Ta).'*flip(d).', (Tc).', (Tat).'*(Tatt).']);
  
    % Compose MEP:
    system = systemstruct(p);
    eigV = [flip(b), g];

    idx = zeros(length(eigV),1);
    for i = 1:length(eigV)
        idx(i) = find(system.map == eigV(i));
    end
    mep = conversion(system,idx);
end

