%% logistic_model.m
%
% This function computes the RHS vector of the logistic function.
%
% Args:
%
% * |x| - vector of state variables [mac mHCO3 mCH4 X]'
% * |p| - vector of parameter values
%
% ret:
% * |rhs| - times at which solution is output

function [rhs] = logistic_RHS(t,x,p)
    
    % specify system of equations
    
    rhs = p(1)*x - p(2)*x^2;

end