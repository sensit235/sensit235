%% scale_parameters.m
%
% A function to scale the experiments generated in the Morris method to the
% parameter values used in the logistic function.

function [p_]=scale_parameters(p_min,p_max,x)

    p_=p_min;
    p_ = (p_max-p_min).*x + p_min;
    
end