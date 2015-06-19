%% tissue_RHS.m
%
% This is an unpublished model of tissue shrinkage as a result of elevated
% temperatures.

function [rhs] = tissue_RHS(t,x,p,T)

k =@(t) exp(p(1)+p(2)/T(t));
l =@(t) exp(p(3)+p(4)/T(t));
m =@(t) exp(p(5)+p(6)/T(t));

rhs = zeros(3,1);
rhs(1) = -k(t)*x(1)+m(t)*x(2);
rhs(2) = -l(t)*x(2)-m(t)*x(2)+k(t)*x(1);
rhs(3) = l(t)*x(2);

end