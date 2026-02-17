% This code reproduces the numerical results reported in Example 5.2 of 
% the accompanying manuscript.

% It uses the MacaulayLab package for constructing and solving 
% multiparameter eigenvalue problems (MEPs). MacaulayLab is available 
% for download at: https://www.macaulaylab.net

%% Initialize problem setup:

clc; clear all; clf; 

% Construct H:
c = flip([1; -1.2024; 2.3675; -2.0039; 2.2337; -1.0420; 0.8513]);
d = flip([0.0448; 0.2368; 0.0013; 0.0211; 0.2250; 0.0219]);
N = length(c)-1;
n = 2;
isContinuous = false;
G = @(z) polyval(flip(d),z)./polyval(flip(c),z);

% Compute H2-norm of H:
fun = @(w) (G(exp(1i*w)).*conj(G(exp(1i*w))));
G_norm = sqrt(1/(2*pi)*integral(@(t) fun(t),-pi,pi))

%% Compute residuals over grid of feasible (a0, a1)-values

% Create grid:
NbPoints = 50;
a0 = linspace(-5, 5, NbPoints);
a1 = linspace(-5, 5, NbPoints);
[X,Y] = meshgrid(a1,a0);

% Iterate over feasible values:
Z = zeros(length(a0), length(a1));
for i = 1:length(a0)
    for j = 1:length(a1)
        a = [a0(i);a1(j);1];

        % Only use stable models:
        if max(abs(roots(flip(a)))) > 1
            Z(i,j) = NaN;
            continue
        end

        % Compute optimal b:
        b = getOptimalZeros(c,d,a,isContinuous);

        % Compute cost function value:
        Gr = @(z) polyval(flip(b),z)./polyval(flip(a),z);
        fun = @(w) (Gr(exp(1i*w))-G(exp(1i*w))).*conj(Gr(exp(1i*w))-G(exp(1i*w))); 
        Z(i,j) = sqrt(1/(2*pi)*integral(@(t) fun(t),-pi,pi));
    end
end

%% Plot contours:
figure(1);
levels = [0.86, 0.88, 0.9, 0.94, 0.98, 1, 1.04,...
    1.08, 1.1, 1.12, 1.16, 1.2, 1.24, 1.26, 1.27, 1.28, 1.285, 1.289, 1.291, 1.2915, 1.2918];
contour(X,Y,Z,levels)
colormap jet
colorbar

%% Construct and solve MEP:

sys = tf(d.',c.',-1); 
[mep, p, M] = getMep(sys,n)
[sols,output] = macaulaylab(mep,30,posdim=true,verbose=true);
sols.values

%% Evaluate Theorem 3.6:
nb_sols = 0;
for j = 0:n
    nb_sols = nb_sols + nchoosek(N-n-1+j,j)*2^j;
end
disp(['N,n: ',mat2str(N),',',mat2str(n),' --> ',mat2str(nb_sols)])
