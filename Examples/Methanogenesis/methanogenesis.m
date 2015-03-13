%% methanogenesis
%
% This function solves the model of methanogenesis given to us by Chewy.
%
% args:
%   t - equivalent of tspan in ode15s - time range of solution
%   p - data structure defined in parameters.m with parameter values
%   i - vector of initial conditions [mac mHCO3 mCH4 X]'
%
% ret:
%   t - times at which solution is output
%   x - solution vector x(# state variable,time)

function [t x] = methanogenesis(t,p,i)
    
    % specify auxiliary equations
    r = @(x) p.k*x(4)*x(1)/(p.Kac + x(1))*...
        (1 - exp(...
        (p.DG0+p.R*p.T*log(x(2)*x(3)/x(1))+p.nup*p.DGp)/...
        p.chi/p.R/p.T));
    
    % specify system of equations
%     s = @(t,x) [-r(x);...
%         r(x);...
%         r(x);...
%         p.Y*r(x) - p.m*x(4)];
    
    s = @(t,x) [-r(x);...
        r(x);...
        r(x);...
        p.Y*r(x)*x(3)/(x(3)+p.Kn) - p.m*x(4)];
    
    % solve system of equations
    options=odeset('RelTol',1e-6,'AbsTol',1e-6,'NormControl','on',...
        'InitialStep',1e-4);
    [t x] = ode15s(s,t,i,options);
%     [t x] = ode45(s,t,i);

end