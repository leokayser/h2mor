% This code reproduces the numerical results reported in Section 5 of 
% the accompanying manuscript.

% It uses the MacaulayLab package for constructing and solving 
% multiparameter eigenvalue problems (MEPs). MacaulayLab is available 
% for download at: https://www.macaulaylab.net

%% Initialize problem setup:

clc; clear all;

% Given model:
c = poly([-2, -2, -3]);
d = poly([4, 5]);
N = length(c)-1;
sys = tf(d,c); 

%% Get and solve MEP:

n = 1;
[mep, p, M] = getMep(sys,n);
[sols,output] = macaulaylab(mep,30,posdim=true,verbose=true);
sols.values

%% Investigate bad locus:

B = M(:, 2:end);

% Trick to reform symbolic matrix polynomial into MEP format:
syms q [N-1 1]
p_bad = B*[1; q]
system_bad = systemstruct(p_bad);
eigV = [q]; % (Eliminate q)
idx = zeros(length(eigV),1);
for i = 1:length(eigV)
    idx(i) = find(system_bad.map == eigV(i));
end
mep_bad = conversion(system_bad,idx);

% Solve mep for bad locus:
[sols,output] = macaulaylab(mep_bad,30,posdim=true,verbose=true);
sols.values

%% Evaluate the solution count:
nb_sols = 0;
for j = 0:n
    nb_sols = nb_sols + nchoosek(N-n-1+j,j)*2^j;
end
nb_sols = nb_sols - nchoosek(N-2, n-1); % Correction term
disp(['N,n: ',mat2str(N),',',mat2str(n),' --> ',mat2str(nb_sols)])