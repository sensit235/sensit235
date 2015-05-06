%% test_function_3.m
%
% This is the Sobol g-function as presented on page 39 of the Saltelli book
% variance [0.7165 0.1791 0.0237 0.0072 0.0001 0.0001 0.0001 0.0001]
%

function [y] = test_function_3(x)

w = x;
g = @(x,a) (abs(4*x-2) + a)/(1+a);
a = [0 1 4.5 9 99 99 99 99];
y = 1;
for i=1:8
    y = y.*g(w(i,:),a(i));
end

end