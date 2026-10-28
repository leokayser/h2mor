% This code reproduces the numerical results reported in Section 5 of 
% the accompanying manuscript.

% It uses the MacaulayLab package for constructing and solving 
% multiparameter eigenvalue problems (MEPs). MacaulayLab is available 
% for download at: https://www.macaulaylab.net

clc; clear all;

%% State-space model:
A = [-1, 1, 0 -3/8; 
    0, -1, 1, 3/4; 
    0, 0, -1, -1/2;
    -1/2, 3/4, -3/8, -10];
B = [0; 0; 1; 1/16];
C = [1, 0, 0, 1/16];

N = size(A,1);
sys = ss(A,B,C,[]);  % continuous time (Ts = 0)
sys = tf(sys)
c = sys.denom{1};
d = sys.num{1}(2:length(d));

%% Get and solve MEP:
n = 3;
[mep, p, M] = getMep(sys,n);

[sols,output] = macaulaylab(mep,30,posdim=true,verbose=true);
sols.values

%% Evaluate the solution count:
nb_sols = 0;
for j = 0:n
    nb_sols = nb_sols + nchoosek(N-n-1+j,j)*2^j;
end
disp(['N,n: ',mat2str(N),',',mat2str(n),' --> ',mat2str(nb_sols)])

%% Filter real and stable models:

realPoints = sols.values(all(abs(imag(sols.values)) <= 1e-6,2),:);

isStable = false(size(realPoints,1),1);
for i = 1:size(realPoints,1)
    isStable(i) = all(real(roots([1, realPoints(i,:)])) < 0);
end
stablePoints = realPoints(isStable, :) % (s+1)^3

%% Solve restricted problem:
% Compose Tp coefficient matrix:
syms p 
a = [-p^3,3*p^2, -3*p,1]; % a(s) = (s-p)^3;
Tp3 = getToeplitz(a, N);

% Compose Tc coefficient matrix:
Tc = getToeplitz(flip(c), 3)

% Compose alternative Ta_tilde matrix:  
tmp = [p^4, 4*p^3, 6*p^2, 4*p, 1] % a_tilde(z)' = (s+p)^4;
Tp4 = getToeplitz(tmp, N-1)

syms b [1 3]
syms g [1 n]

% Get system of equations:
p = (Tp3).'*flip(d).' - (Tc).'*flip(b).' - (Tp4).'*g.';
M = vpa([(Tp3).'*flip(d).', (Tc).', (Tp4).']);

% Compose MEVP:
system = systemstruct(p);
eigV = [flip(b), g];

idx = zeros(length(eigV),1);
for i = 1:length(eigV)
    idx(i) = find(system.map == eigV(i));
end
mep = conversion(system,idx);

dend = 50;
[sols,output] = macaulaylab(mep,dend,posdim=true,verbose=true);
critPoints = sols.values; % Optimal values of pole p

% Filter real solutions only:
realPoints = critPoints(abs(imag(critPoints)) <= 1e-6)

% Stable solutions:
stabPoints = realPoints(realPoints < 0)
