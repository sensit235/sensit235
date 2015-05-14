%% methanogenesis.m
%
% Analysis of Chewy's methanogenesis model.

close all
clear all
addpath ../../Methods/Morris' Method'/
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

morris_methanogenesis_RHS = @(x,theta0)...
    methanogenesis_RHS(x,...
    [DG0, R, theta0(1), theta0(2), DGp, theta0(3), theta0(4), T,...
    theta0(5), Kn, theta0(6)]);

%%
% For the global sensitivity analysis the parameters will be varied by a
% fixed percentage |pcg| of their reference value.
%
% * theta_min = vector of lower bounds of parameters under consideration
% * theta_max = vector of upper bounds of parameters under consideration

theta = [k, nup, chi, Y, Kac, m];
pcg = 0.0001;
theta_min = (1 + pcg/100).*theta;
theta_max = (1 - pcg/100).*theta;

%%
% The Morris method can now be used to compute the sensitivity measures for
% the model parameters of interest.

tspan = linspace(0,20*24*60*60,100);
[mnt0 sdt0] = ...
    run_morris(4,4,tspan,x0,theta_min,theta_max,morris_methanogenesis_RHS);

%%
% Now results are obtained without normalisation and with the same
% normalisation as the relative sensitivity functions. And are plotted
% against the analytic solutions.

[mnt1 sdt1] = ...
    run_morris(4,4,tspan,x0,theta_min,theta_max,morris_methanogenesis_RHS,'none');

[mnt2 sdt2] = ...
    run_morris(4,4,tspan,x0,theta_min,theta_max,morris_methanogenesis_RHS,'rsf');

leg_text = {'k', 'nup', 'chi', 'Y', 'Kac', 'm'};

% standard morris
figure
plot(tspan,mnt0')
% hold on
% plot(tspan,mnt1')
title('Mean vs Time','Interpreter','LaTex','FontSize',20)
set(gca,'FontSize',14)

legend(leg_text,'Interpreter','LaTex','FontSize',20)

print -depsc figure1

% no normalisation
figure
% plot(tspan,mnt0')
% hold on
plot(tspan,mnt1')
title('Mean vs Time','Interpreter','LaTex','FontSize',20)
set(gca,'FontSize',14)

legend(leg_text,'Interpreter','LaTex','FontSize',20)

print -depsc figure2

% rsf normalisation
figure
% plot(tspan,mnt0')
% hold on
plot(tspan,mnt2')
title('Mean vs Time','Interpreter','LaTex','FontSize',20)
set(gca,'FontSize',14)

legend(leg_text,'Interpreter','LaTex','FontSize',20)

print -depsc figure3

%% Computing the sensitivity measures using total sensitivity functions.
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

print -depsc figure4

%% Comparison of the two methods without normalisation.

figure
plot(tspan,mnt1',t,y(:,5:10),'blackx')
title('Mean vs Time','Interpreter','LaTex','FontSize',20)
set(gca,'FontSize',14)

legend(leg_text,'Interpreter','LaTex','FontSize',20)

print -depsc figure5