%% test_jacobian.m
%
% Script exploring the DERIVESTsuite to get numerical jacobian.

addpath ../../External/DERIVESTsuite


%%
% Specify the parameters

DG0 = -15802.1961; % standard value of Gibbs free energy at T=310.15K (J/mol)
R = 8.3145; % gas constant (J/mol/K)
k = 2.5e-6; % (mol/g/s) **
nup = 0.5; % (per reaction) **
DGp = 45000; % phosphorylation energy (J/(mol ATP))
chi = 2; % (per reaction) **
Y = 2.1; % (Ymax) (g/mol) **
T = 310.15; % physiological temperature (K)
Kac = 5e-3; % (Kd) half saturation constant (molal) **
Kn = 0; % effect of nutrient (molal)
m = 2.2e-7; % (D) specific maintenance rate (1/s) **

%%
% Define the system of equations for Jacobian

% specify auxiliary equations
r = @(x) k*x(4)*x(1)/(Kac + x(1))*...
    (1 - exp(...
    (DG0+R*T*log(x(2)*x(3)/x(1))+nup*DGp)/...
    chi/R/T));

% specify system
s = @(t,x) [-r(x);...
    r(x);...
    r(x);...
    Y*r(x)*x(3)/(x(3)+Kn) - m*x(4)];

x0 = [3.5e-3, 45.2e-3, 1e-6, 0.009];
theta = [k, nup, chi, Y, Kac, m];

model = @(x) s(0,x);
[grad,err] = jacobianest(model,x0)

%%
% Define the system of equations for derivative wrt theta

% specify auxiliary equations
r = @(x,p) p(1)*x(4)*x(1)/(p(5) + x(1))*...
    (1 - exp(...
    (DG0+R*T*log(x(2)*x(3)/x(1))+p(2)*DGp)/...
    p(3)/R/T));

% specify system
s = @(t,x,p) [-r(x,p);...
    r(x,p);...
    r(x,p);...
    p(4)*r(x,p)*x(3)/(x(3)+Kn) - p(6)*x(4)];

x0 = [3.5e-3, 45.2e-3, 1e-6, 0.009];
theta = [k, nup, chi, Y, Kac, m];

model = @(x) s(0,x0,x);
[grad,err] = jacobianest(model,theta)