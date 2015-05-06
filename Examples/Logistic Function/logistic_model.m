%% logistic
%
% This function solves the logistic function.
%
% args:
%   t - equivalent of tspan in ode15s - time range of solution
%   p - vector of parameters
%   i - vector of initial conditions
%
% ret:
%   t - times at which solution is output
%   x - solution vector x(# state variable,time)

function [t x] = logistic(t,p,i)
    
    % specify system of equations
    
    s = @(t,x) p(1)*x - p(2)*x^2;
    
    % solve system of equations
    options=odeset('RelTol',1e-6,'AbsTol',1e-6,'NormControl','on',...
        'InitialStep',1e-4);
    [t x] = ode15s(s,t,i,options);

end