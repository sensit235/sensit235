function [t_num,y_num] = sensit_tsf(t,model_rhs,theta,x0)
%% run_tsf.m
%
% This function computes the total sensitivity functions for a given model
% by numerically solving the sensitivity equations.
%
% Args:
%
% * |t| - vector of times at which to compute sensitivity
% * |model| - a handle to the model
% * |theta| - vector of optimal parameter values
% * |x0| - vector of initial conditions


addpath ../../External/DERIVESTsuite/

n_p = length(theta); % number of parameters
n_s = length(x0); % number of states

%%
% Compute the solution to the sensitivity equations constructed in
% tsf_model.

xinit = zeros(n_s + n_s*n_p +n_s*n_s,1);
xinit(1:n_s) = x0;
xinit(n_s+1:n_s+n_s*n_p) = zeros(n_s*n_p,1); % add initial condition for param sens
xic = eye(n_s); % initial conditions for initial condition sens
xinit(n_s+n_s*n_p+1:n_s + n_s*n_p +n_s*n_s) = reshape(xic,n_s*n_s,1); % full initial condition vector
options = odeset('RelTol',1e-6,'AbsTol',1e-6);
[t_num,y_num] = ode15s(@tsf_model,t,xinit,options,model_rhs,theta,n_s,n_p);

end

%% tsf_model
%
% This function returns the rhs required to compute the total sensitivity
% functions and requires the rhs of the model as input.
%
% Args:
%
% * |t| - time
% * |v| - vector representing RHS of sensitivity functions
% * |model| - function handle to RHS of model of interest
% * |pars| - model parameters
% * |n_s| - number of states
% * |n_p| - number of parameters

function dy = tsf_model(t,v,model,pars,n_s,n_p)

%%
% Compute the current value of RHS for state variables

dy_state = model(t,v(1:n_s),pars);

%%
% compute the Jacobian df/dx

dfdx = jacobianest(@(x) model(t,x,pars),v(1:n_s));

%%
% compute the derivate $\frac{dF}{d\theta}$

dfda = jacobianest(@(x) model(t,v(1:n_s),x),pars);

%%
% compute the sensitivity matrix
    
v40     = v(n_s+1:n_s+n_s*n_p);      % the sensitiity coordinates
vv      = reshape(v40,[n_p,n_s])';   % vector 2 matrix 
SenEq   = dfdx*vv+dfda;              % SENSITIVITY EQUATION
dy_sens = reshape(SenEq',n_s*n_p,1); % matrix 2 vector

%%
% add initial conditions

v_IC       = v(n_s+n_s*n_p+1:n_s+n_s*n_p+n_s*n_s);
vv_IC      = reshape(v_IC,[n_s,n_s])';
Sen_IC     = dfdx*vv_IC;
dy_sens_IC = reshape(Sen_IC',n_s*n_s,1);

%%
% append sensitivity equations as vector

dy      = [dy_state; dy_sens; dy_sens_IC];            % returning full eqns

end