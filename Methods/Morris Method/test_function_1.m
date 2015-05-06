%% test_function_1.m
%
% This is the test function described in equation 10 of the paper "The use
% of graph theory in the sensitivity analysis of the model output: a
% second order screening method", Campolongo & Braddock 1999
%
% This paper has been corrected in 2002, and analysis should be compared to
% that publication.
%
% first order affect [90.050 71.045 41.470 20.825]
% standard deviation [31.08 28.86 17.75 11.54]

function y = test_function_1(x)

% w = (x+1)/2;
w = x;
bi = [0.05 0.59 10 0.21];
bij = [0 80 60 40; 0 30 0.73 0.18; 0 0 0.64 0.93; 0 0 0 0.06];
y = bi*w + diag(w'*bij*w)';

end