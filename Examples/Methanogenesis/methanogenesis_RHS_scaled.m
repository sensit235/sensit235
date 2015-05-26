%% methanogenesis_RHS.m
%
% This function computes the RHS vector of the model of methanogenesis
% given to us by Chewy.
%
% Args:
%
% * |x| - vector of state variables [mac mHCO3 mCH4 X]'
% * |p| - vector of parameter values
%
% ret:
% * |rhs| - times at which solution is output

function [rhs] = methanogenesis_RHS_scaled(x,p)
    
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
        (DG0+R*T*log(x(2)*x(3)/1000/x(1))+nup*DGp)/...
        chi/R/T));
    
    rhs = [-r(x);...
        r(x);...
        r(x)*1000;...
        Y*r(x)*x(3)/1000/(x(3)/1000+Kn) - m*x(4)];

end