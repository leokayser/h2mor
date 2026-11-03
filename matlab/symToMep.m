function mep = symToMep(eqs, linVars)
%SYMTOMEP  MacaulayLab MEP from polynomial equations that are linear in some variables.
%   mep = symToMep(eqs, linVars) converts the symbolic equations eqs, which
%   are linear in the variables linVars, into a MacaulayLab MEP. The
%   remaining variables of eqs become the eigenvalues.

    system = systemstruct(eqs);
    idx = zeros(numel(linVars), 1);
    for i = 1:numel(linVars)
        idx(i) = find(system.map == linVars(i));
    end
    mep = conversion(system, idx);
end