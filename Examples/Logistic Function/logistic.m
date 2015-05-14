
%% logistic.m
%
% The logistic equation detailed in |\cite{Banks2007a}| provides a good
% test of the performance of the Morris method. It has an analytic solution
% and few parameters allowing the direct computation of errors and
% interpretation of results.

close all
clear all
arun_sens_logistic_explicit % generate analytic solution
addpath ../../Methods/Morris' Method'/
addpath ../../Methods/TSF/

%% Define the parameters and initial conditions.
%
% * theta = optimal fits of parameters under consideration $\theta_0$
% * x0 = initial conditions under consideration

theta = [0.8,0.1];
x0 = [0.3];
tspan = [0 16];

%% Computing the sensitivity measures using the Morris method.
%
% For the global sensitivity analysis the parameters will be varied by a
% fixed percentage |pcg| of their reference value.
%
% * theta_min = vector of lower bounds of parameters under consideration
% * theta_max = vector of upper bounds of parameters under consideration

pcg = 0.0001;
theta_min = (1 + pcg/100).*theta;
theta_max = (1 - pcg/100).*theta;
x0_min = (1 + pcg/100).*x0;
x0_max = (1 - pcg/100).*x0;

%%
% The Morris method can now be used to compute the sensitivity measures for
% the model parameters.

tspan = linspace(tspan(1),tspan(2),100);
[mnt0 sdt0] = ...
    run_morris(4,4,tspan,x0_min,x0_max,theta_min,theta_max,@logistic_RHS);

%%
% Now results are obtained without normalisation and with the same
% normalisation as the relative sensitivity functions. And are plotted
% against the analytic solutions.

[mnt1 sdt1] = ...
    run_morris(4,4,tspan,x0_min,x0_max,theta_min,theta_max,@logistic_RHS,'none');

[mnt2 sdt2] = ...
    run_morris(4,4,tspan,x0_min,x0_max,theta_min,theta_max,@logistic_RHS,'rsf');

% generate analytic solutions
load ref t dx_da dx_db dx_dc Ndx_da Ndx_db Ndx_dc
leg_text = {'a','b','x0'};

% standard morris vs no normalisation and analytic
figure
plot(tspan,mnt1')
hold on
plot(tspan,mnt0')
title('Mean vs Time','Interpreter','LaTex','FontSize',20)
set(gca,'FontSize',14)
plot(t,dx_da,'blackx')
plot(t,dx_db,'blackx')
plot(t,dx_dc,'blackx')
legend(leg_text,'Interpreter','LaTex','FontSize',20)

print -depsc figure1

% standard morris vs rsf normalisation and analytic
figure
plot(tspan,mnt2')
hold on
plot(tspan,mnt0')
title('Mean vs Time','Interpreter','LaTex','FontSize',20)
set(gca,'FontSize',14)
plot(t,Ndx_da,'blackx')
plot(t,Ndx_db,'blackx')
plot(t,Ndx_dc,'blackx')
legend(leg_text,'Interpreter','LaTex','FontSize',20)

print -depsc figure2

%% Compute the sensitivity measures using total sensitivity functions.

%%
% The vector of parameters under consideration is constructed and the
% |run_tsf| function called to compute the sensitivities

[tl yl] = run_tsf(tspan,@logistic_RHS,theta,x0);

% standard tsf vs analytic
figure
plot(tl,yl(:,2:3))
hold on
title('Mean vs Time','Interpreter','LaTex','FontSize',20)
set(gca,'FontSize',14)
plot(t,dx_da,'blackx')
plot(t,dx_db,'blackx')
plot(t,dx_dc,'blackx')
legend(leg_text,'Interpreter','LaTex','FontSize',20)

print -depsc figure3