% This code reproduces the numerical results reported in Section 5 of 
% the accompanying manuscript.

% It uses the MacaulayLab package for constructing and solving 
% multiparameter eigenvalue problems (MEPs). MacaulayLab is available 
% for download at: https://www.macaulaylab.net

clc; clear all;

%% Given model:
% Paper: Globally Optimal H2-Norm Model Reduction: A Numerical
% Linear Algebra Approach, O.M. Agudelo, V. Christof, B. De Moor
% (2021) - IFAC-PapersOnLine.
d = [1, 9, -10];        % numerator
c = [1, 12, 49, 78];    % denominator, poles at -6 and -3 +/- 2i
N = length(c)-1;
sys = tf(d,c);      % continuous time (Ts = 0)

%% Get and solve MEP:
n = 1;
[mep, p, M] = getMep(sys,n);
roots(flip(coeffs(det(M))))
[sols,output] = macaulaylab(mep,30,posdim=true,verbose=true);
sols.values

%% Evaluate the solution count:
nb_sols = 0;
for j = 0:n
    nb_sols = nb_sols + nchoosek(N-n-1+j,j)*2^j;
end
disp(['N,n: ',mat2str(N),',',mat2str(n),' --> ',mat2str(nb_sols)])

%% Plot the optimization landscape:
% Create grid:
NbPoints = 100;
a0 = linspace(0.05, 12, NbPoints);
H = @(s) polyval(d,s)./polyval(c,s);

% Iterate over feasible values:
costs = zeros(1,length(a0));
for i = 1:length(a0)
    a = [1,a0(i)];
    b = getOptimalZeros(flip(c).',flip(d).',flip(a).',true);
    Hr = @(s) polyval(b,s)./polyval(a,s);
    fun = @(w) (H(1i*w)-Hr(1i*w)).*conj(H(1i*w)-Hr(1i*w)); 
    costs(i) = sqrt(1/(2*pi)*integral(@(t) fun(t),-inf,inf));
end

plot(a0,costs)