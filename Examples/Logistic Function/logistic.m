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

%% Define the parameters and initial conditions.
%
% * theta = optimal fits of parameters under consideration $\theta_0$
% * x0 = initial conditions under consideration

theta = [0.8,0.1];
x0 = [0.3];

%% Implementing the model.
%
% The actual ode / system of odes must be implemented as a function, in
% this case |logistic_model.m|. The model function must take the parameters
% , initial conditions and a time span as the arguments and return a
% solution _y_ to the model.

tspan = [0 16];
[t x] = logistic_model(tspan,theta,x0);
plot(t,x)
snapnow;

%% Computing the sensitivity measures using the Morris method.
%
% The Morris method has no concept of initial conditions being distinct
% from parameters and therefore we need to define a new function that can
% deal with this.

morris_logistic_model = @(t,theta)...
    logistic_model(t,theta(1:2),theta(3));

%%
% A new parameter vector is now defined

theta = [theta x0];

%%
% For the global sensitivity analysis the parameters will be varied by a
% fixed percentage |pcg| of their reference value.
%
% * theta_min = vector of lower bounds of parameters under consideration
% * theta_max = vector of upper bounds of parameters under consideration

pcg = 10;
theta_min = (1 + pcg/100).*theta;
theta_max = (1 - pcg/100).*theta;

%%
% The Morris method can now be used to compute the sensitivity measures for
% the model parameters.

tspan = linspace(tspan(1),tspan(2),100);
[mnt0 sdt0] = ...
    run_morris(4,4,tspan,theta_min,theta_max,morris_logistic_model);

%%
% Now results are obtained without normalisation and with the same
% normalisation as the relative sensitivity functions. And are plotted
% against the analytic solutions.

[mnt1 sdt1] = ...
    run_morris(4,4,tspan,theta_min,theta_max,morris_logistic_model,'none');

[mnt2 sdt2] = ...
    run_morris(4,4,tspan,theta_min,theta_max,morris_logistic_model,'rsf');

% generate analytic solutions
load ref t dx_da dx_db dx_dc Ndx_da Ndx_db Ndx_dc
leg_text = {'a','b','x0'};

% standard morris vs no normalisation and analytic
figure
plot(tspan,mnt0')
hold on
plot(tspan,mnt1')
title('Mean vs Time','Interpreter','LaTex','FontSize',20)
set(gca,'FontSize',14)
plot(t,dx_da,'blackx')
plot(t,dx_db,'blackx')
plot(t,dx_dc,'blackx')
legend(leg_text,'Interpreter','LaTex','FontSize',20)

print -depsc figure1

% standard morris vs rsf normalisation and analytic
figure
plot(tspan,mnt0')
hold on
plot(tspan,mnt2')
title('Mean vs Time','Interpreter','LaTex','FontSize',20)
set(gca,'FontSize',14)
plot(t,Ndx_da,'blackx')
plot(t,Ndx_db,'blackx')
plot(t,Ndx_dc,'blackx')
legend(leg_text,'Interpreter','LaTex','FontSize',20)

print -depsc figure2