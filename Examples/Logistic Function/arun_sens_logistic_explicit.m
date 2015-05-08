clear all
close all

% Parameters 
% 0.8 0.01 0.1
a  = 0.8;
b  = 0.1;
x0 = 0.3;   
k = (a/x0)-b;

% Interval
t_end = 16;
t = 0:1:t_end;

% Solution
x_sol = a./(b + (-b + a/x0)./exp(a.*t));

dx_da = (exp(a.*t).*x0.*(b*(-1+exp(a.*t))*x0 + a.*t.*(a - b.*x0)))./(a + b*(-1 + exp(a.*t)).*x0).^2;
dx_db = -((a.*exp(a.*t).*(-1 + exp(a.*t)).*x0.^2)./(a+b.*(-1 + exp(a.*t)).*x0).^2);

dx_dc = a^2*exp(a.*t)./(a+b*(exp(a.*t)-1)*x0).^2;

% Normalization
Ndx_da =  dx_da.*(a./x_sol);
Ndx_db =  dx_db.*(b./x_sol);
Ndx_dc =  dx_dc.*(x0./x_sol);

save ref.mat