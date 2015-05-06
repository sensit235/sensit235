%% logistic.m
%
% The logistic equation detailed in |\cite{Banks2007a}| provides a good
% test of the performance of the Morris method. It has an analytic solution
% and few parameters allowing the direct computation of errors and
% interpretation of results.

%% Define the parameters and initial conditions.
%
% * theta = optimal fits of parameters under consideration $\theta_0$
% * theta_min = vector of lower bounds of parameters under consideration
% * theta_max = vector of upper bounds of parameters under consideration
% * x0 = initial conditions under consideration

theta = [0.8,0.1];
x0 = [0.3];

% In this case we will vary the parameters by a fixed percentage |pcg| of
% their mean value to explore the Morris method.

pcg = 10;
theta_min = (1 + pcg/100).*theta;
theta_max = (1 - pcg/100).*theta;

%% Implementing the model.
%
% The actual ode / system of odes must be implemented as a function, in
% this case |logistic_model.m|. The model function must take the parameters
% , initial conditions and a time span as the arguments and return a
% solution _y_ to the model.