%% methanogenesis_tsf.m
%
% This script shows how to compute the total sensitivity functions for the
% methanogenesis example using run_tsf.

close all
clear all
addpath ../../Methods/TSF/

%% Define the parameters and initial conditions.
%
% * theta = optimal fits of parameters under consideration $\theta_0$
% * x0 = initial conditions under consideration

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

% parameter vector
theta = [DG0, R, k, nup, DGp, chi, Y, T, Kac, Kn, m];

% initial condition
x0 = [3.5e-3, 45.2e-3, 1e-6, 0.009]; % [molal molal molal g/kg]

%% Computing the sensitivity measures using the Morris method.
%
% Need to define a wrapper function that only uses a subset of the actual
% parameter values in the analysis.

tsf_methanogenesis_RHS = @(x,theta)...
    methanogenesis_RHS(x,...
    [DG0, R, theta(1), theta(2), DGp, theta(3), theta(4), T,...
    theta(5), Kn, theta(6)]);

%% Run the sensitivity analysis.
%
% The vector of parameters under consideration is constructed and the
% |run_tsf| function called to compute the sensitivities

theta = [k, nup, chi, Y, Kac, m];
tspan = linspace(0,20*24*60*60,100);
[t y] = run_tsf(tspan,tsf_methanogenesis_RHS,theta,x0);

%%
% Plot the results for mac

leg_text = {'k', 'nup', 'chi', 'Y', 'Kac', 'm'};
figure
plot(t,y(:,5:10))
title('Mean vs Time','Interpreter','LaTex','FontSize',20)
set(gca,'FontSize',14)

legend(leg_text,'Interpreter','LaTex','FontSize',20)

print -depsc figure2