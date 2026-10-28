% This code reproduces the numerical results reported in Section 5 of 
% the accompanying manuscript.

% It uses the MacaulayLab package for constructing and solving 
% multiparameter eigenvalue problems (MEPs). MacaulayLab is available 
% for download at: https://www.macaulaylab.net

clc; clear all;

%% MOR contour plot

% Four disk system:
d = [0.0448, 0.2368, 0.0013, 0.0211, 0.2250, 0.0219];
c = [1, -1.2024, 2.3675, -2.0039, 2.2337, -1.0420, 0.8513];
fourDisk = tf(d,c,-1)
N = length(c)-1;
isContinuous = false;

% Compose Tp coefficient matrix:
syms p 
a = [p^2,-2*p,1]; % a(z) = (z-p)^2;
Tp2 = getToeplitz(a, N);

% Compose Tc coefficient matrix:
Tc = getToeplitz(flip(c), 2)

% Compose alternative Ta_tilde matrix:  
tmp = flip([-p^3, 3*p^2, -3*p, 1]) % a_tilde(z) = (-pz + 1)^3;
Tp3 = getToeplitz(tmp, N-1)

syms b [1 2]
syms g [1 n-1]

% Get system of equations:
p = (Tp2).'*flip(d).' - (Tc).'*flip(b).' - (Tp3).'*g.';
M = vpa([(Tp2).'*flip(d).', (Tc).', (Tp3).']);

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
realPoints = critPoints(abs(imag(critPoints)) <= 1e-6); 

% For each critical point (value of the pole) find a(z):
models = zeros(length(realPoints),3);
for i = 1:length(realPoints)
    models(i,:) = flip(poly([realPoints(i), realPoints(i)])); % [a0, a1, a2] with a0 + a1z + a2z^2
end

% Filter solutions based on the region displayed in the figure:
models_disp = models(abs(models(:,2)) < 2 & abs(models(:,1)) < 1,:)

%% Solve the n=1 problem:
sys = tf(flip(d),flip(c),-1); 
[mep, p, M] = getMep(sys,1)
[sols,output] = macaulaylab(mep,30,posdim=true,verbose=true);
critPoints = sols.values
poles = -1./critPoints
% filter real and stable models
poles = poles(abs(imag(poles)) < 1e-6);
poles = poles(abs(poles)<1)
% These are the local maximisers on the discriminant