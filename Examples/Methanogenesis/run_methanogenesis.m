%% run_methanogenesis.m
%
% Analysis of Chewy's methanogenesis model.

close all
clear all
addpath ../../Methods/Morris' Method'/
addpath ../../Methods/TSF/

%% Define the parameters and initial conditions.

DG0 = -15802.1961;  % Gibbs free energy at T=310.15K (J/mol)
R   = 8.3145;       % gas constant (J/mol/K)
k   = 2.5e-6;       % (mol/g/s) **
nup = 0.5;          % (per reaction) **
DGp = 45000;        % phosphorylation energy (J/(mol ATP))
chi = 2;            % (per reaction) **
Y   = 2.1;          % (Ymax) (g/mol) **
T   = 310.15;       % physiological temperature (K)
Kac = 5e-3;         % (Kd) half saturation constant (molal) **
Kn  = 0;            % effect of nutrient (molal)
m   = 2.2e-7;       % (D) specific maintenance rate (1/s) **

% parameter vector
theta = [DG0, R, k, nup, DGp, chi, Y, T, Kac, Kn, m];

% initial condition
x0 = [3.5e-3, 45.2e-3, 1e-6, 0.009]; % [molal molal molal g/kg]

% A wrapper function that only uses a subset of the actual
% parameter values in the analysis.

restrict_methanogenesis_RHS = @(t,x,theta0)...
    methanogenesis_RHS(x,...
    [DG0, R, theta0(1), theta0(2), DGp, theta0(3), theta0(4), T,...
    theta0(5), Kn, theta0(6)]);

% For the global sensitivity analysis the parameters will be varied by a
% fixed percentage |pcg| of their reference value.
%
% * theta_min = vector of lower bounds of parameters under consideration
% * theta_max = vector of upper bounds of parameters under consideration

theta       = [k, nup, chi, Y, Kac, m];
pcg         = 1e-1;
theta_min   = (1 + pcg/100).*theta;
theta_max   = (1 - pcg/100).*theta;
x0_min      = (1 + pcg/100).*x0;
x0_max      = (1 - pcg/100).*x0;

%% Computing the sensitivity measures using the Morris method

% The Morris method (standard normalization) 
tspan           = linspace(0,20*24*60*60,100);
[mnt0, sdt0]    = ...
    sensit_morris(4,4,tspan,x0_min,x0_max,theta_min,theta_max,restrict_methanogenesis_RHS);

% Morris method (no normalisation) 
[mnt1, sdt1] = ...
    sensit_morris(4,4,tspan,x0_min,x0_max,theta_min,theta_max,restrict_methanogenesis_RHS,'none');

% Morris method (normalization as in relative sensitivity functions)
[mnt2, sdt2] = ...
    sensit_morris(4,4,tspan,x0_min,x0_max,theta_min,theta_max,restrict_methanogenesis_RHS,'rsf');


%% Computing the sensitivity measures using Total Sensitivity Functions.
%
% The vector of parameters under consideration is constructed and the
% |sensit_tsf| function called to compute the sensitivities

theta = [k, nup, chi, Y, Kac, m];
tspan = linspace(0,20*24*60*60,100);

[t, y] = sensit_tsf(tspan,restrict_methanogenesis_RHS,theta,x0);


save saved_data.mat mnt0 sdt0 mnt1 sdt1 mnt2 sdt2 t y tspan 
%% Addressing numerical instabilities.
%
% In the previous figure it can be seen that there is numerical instability
% in the computation of the sensitivity of the initial conditions. This
% arises due to the finite precision available in matlab, but can be
% circumvented in many cases.
%
% In this example a simple rescaling of the state variables is sufficient
% to remove these instabilities. If the Morris method results are then
% recomputed using the scaled equations, the sensitivities will be
% identical and the instabilities removed. The steady state values of the
% first 3 state variables have been used for scaling because they can be
% obtained analytically. The fourth state variable has a steady state value
% of 0 and is therefore not scaled to avoid singularities.

% analytic steady-state
ss1=0.0045;
ss2=0.0442;
ss3=0.0010;
ss4=0;

% initial condition
x0 = [(3.5e-3)/ss1, (45.2e-3)/ss2, (1e-6)/ss3, 0.009]; % [molal molal molal g/kg]
x0_min = (1 + pcg/100).*x0;
x0_max = (1 - pcg/100).*x0;

%%
% The Morris method can now be used to compute the sensitivity measures for
% the model parameters of interest.

restrict_methanogenesis_RHS = @(t,x,theta0)...
    methanogenesis_RHS_scaled(x,...
    [DG0, R, theta0(1), theta0(2), DGp, theta0(3), theta0(4), T,...
    theta0(5), Kn, theta0(6)]);
tspan = linspace(0,20*24*60*60,100);
[mnt3, sdt3] = ...
    sensit_morris(4,4,tspan,x0_min,x0_max,theta_min,theta_max,restrict_methanogenesis_RHS,'none');

%%
leg_text = {'k', 'nup', 'chi', 'Y', 'Kac', 'm'};

% Remove the scaling from the solutions
mnt3(11,:) = mnt3(11,:);
mnt3(12,:) = mnt3(12,:)*ss1/ss2;
mnt3(13,:) = mnt3(13,:)*ss1/ss3;
mnt3(14,:) = mnt3(14,:)*ss1;

figure
plot(tspan,mnt3(11:14,:)',t,y(:,29:32),'blackx')
title('Morris vs TSF','Interpreter','LaTex','FontSize',20)
xlabel('Time (s)','Interpreter','LaTex','FontSize',20)
ylabel('$\frac{\partial x}{\partial x_0}$','Interpreter','LaTex','FontSize',20)
set(gca,'FontSize',14)

legend(leg_text,'Interpreter','LaTex','FontSize',20)

%print -depsc figure7