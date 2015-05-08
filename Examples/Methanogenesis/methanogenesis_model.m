%% methanogenesis_model.m
%
% This function solves the model of methanogenesis given to us by Chewy.
%
% Args:
%
% * |t| - equivalent of tspan in ode15s - time range of solution
% * |p| - data structure defined in parameters.m with parameter values
% * |i| - vector of initial conditions [mac mHCO3 mCH4 X]'
%
% ret:
% * |t| - times at which solution is output
% * |x| - solution vector x(# state variable,time)

function [t x] = methanogenesis(t,p,i)
    
    DG0 = p(1); % standard value of Gibbs free energy at T=310.15K (J/mol)
    R = p(2); % gas constant (J/mol/K)
    k = p(3); % (mol/g/s) **
    nup = p(4); % (per reaction) **
    DGp = p(5); % phosphorylation energy (J/(mol ATP))
    chi = p(6); % (per reaction) **
    Y = p(7); % (Ymax) (g/mol) **
    T = p(8); % physiological temperature (K)
    Kac = p(9); % (Kd) half saturation constant (molal) **
    Kn = p(10); % effect of nutrient (molal)
    m = p(11); % (D) specific maintenance rate (1/s) **

    % specify auxiliary equations
    r = @(x) k*x(4)*x(1)/(Kac + x(1))*...
        (1 - exp(...
        (DG0+R*T*log(x(2)*x(3)/x(1))+nup*DGp)/...
        chi/R/T));
    
    s = @(t,x) [-r(x);...
        r(x);...
        r(x);...
        Y*r(x)*x(3)/(x(3)+Kn) - m*x(4)];
    
    % solve system of equations
    options=odeset('RelTol',1e-6,'AbsTol',1e-6,'NormControl','on',...
        'InitialStep',1e-4);
    [t x] = ode15s(s,t,i,options);

end