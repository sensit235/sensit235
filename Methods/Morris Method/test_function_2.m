%% test_function_2.m
%
% This is the test function described in equation 10 of the paper "The use
% of graph theory in the sensitivity analysis of the model output: a
% second order screening method", Campolongo & Braddock 1999
%
% reduced version without matrix, compare to 2002 paper
%
% morris method mean should reproduce bi below
%

function y = test_function_2(x)

% w = (x+1)/2;
w = x;
bi = [0.05 0.59 10 0.21];
y = bi*w;

end