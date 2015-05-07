%% logistic_model.m
%
% This function solves the logistic equation.
%
% Arguments:
%
% * t - equivalent of tspan in ode15s - time range of solution
% * p - vector of parameters
% * i - vector of initial conditions
%
% Returns:
%
% * t - times corresponding to solution vector
% * x - solution vector

function [t x] = logistic_model(t,p,i)
    
    % specify system of equations
    
    s = @(t,x) p(1)*x - p(2)*x^2;
    
    % solve system of equations
    options=odeset('RelTol',1e-6,'AbsTol',1e-6,'NormControl','on',...
        'InitialStep',1e-4);
    [t x] = ode15s(s,t,i,options);

end