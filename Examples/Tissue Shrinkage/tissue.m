%% tissue.m
%
% This script computes the sensitivity of the tissue contraction model.
% This model has an observation function, which has not been implemented
% yet.

close all
clear all
addpath ../../Methods/TSF/
addpath ../../Methods/Morris' Method'/


%% Define the parameters and initial conditions.
%
% * theta = optimal fits of parameters under consideration $\theta_0$
% * x0 = initial conditions under consideration

% compute initial params
theta0=[1.4735    0.0286  0.4306    0.3340];
r_low_k = 1e-6;
bk = log((r_low_k)/theta0(1))/(1/310 - 1/(273+85));
ak = log(r_low_k) - bk/310;
r_low_l = 1e-6;
bl = log((r_low_l)/theta0(2))/(1/310 - 1/(273+85));
al = log(r_low_l) - bl/310;

r_low_m = 0.001;
r_high = 0.5;
bm = log((r_low_m)/r_high)/(1/310 - 1/(273+85));
am = log(r_low_m) - bm/310;

theta = [ak bk al bl am bm theta0(end-1:end)];

% initial condition
x0 = [1 0 0]; % [N U D]

%% Run the tsf sensitivity analysis.
%
% The vector of parameters under consideration is constructed and the
% |run_tsf| function called to compute the sensitivities

T =@(t) (358-(273+37))*heaviside(-t+120)+273+37;
tsf_tissue_RHS =@(t,x,p) tissue_RHS(t,x,p,T);

tspan = [0:1500]';
[t y] = run_tsf(tspan,tsf_tissue_RHS,theta,x0);

%% Compute sensitivities of observation function and plot.
ob1 = [-100,-100*theta(end-1),-100*theta(end)];

% parameters in system of ODEs
S = {};
for i=4:9
    S = [S ob1*(y(:,i:8:27)')];
end

% parameters in observation function
S = [S -100*(y(:,2)')];
S = [S -100*(y(:,3)')];

figure
plot(tspan,100+ob1*([y(:,1) y(:,2) y(:,3)]'))

figure
plot(tspan,[(S{1}')*theta(1) (S{2}')*theta(2) (S{3}')*theta(3)...
    (S{4}')*theta(4) (S{5}')*theta(5) (S{6}')*theta(6)...
    (S{7}')*theta(7) (S{8}')*theta(8)])

%% Run the morris sensitivity analysis.
%
% The vector of parameters under consideration is constructed and the
% |run_tsf| function called to compute the sensitivities

pcg = 1e-1;
theta_min = (1 + pcg/100).*theta;
theta_max = (1 - pcg/100).*theta;
x0_min = (1 + pcg/100).*x0;
x0_max = (1 - pcg/100).*x0;

T =@(t) (358-(273+37))*heaviside(-t+120)+273+37;
morris_tissue_RHS =@(t,x,p) tissue_RHS(t,x,p,T);

tspan = [0:1500]';

[mnt1 sdt1] = ...
    run_morris(4,4,tspan,x0_min,x0_max,theta_min,theta_max,morris_tissue_RHS,'none');

%% Compute sensitivities of observation function and plot.
ob1 = [-100,-100*theta(end-1),-100*theta(end)];

% parameters in system of ODEs
S = {};
for i=4:9
    S = [S ob1*mnt1(i:11:36,:)];
end

% % parameters in observation function
S = [S -100*(mnt1(2,:))];
S = [S -100*(mnt1(3,:))];

figure
plot(tspan,100+ob1*mnt1(1:3,:))

figure
plot(tspan,[(S{1}') (S{2}') (S{3}')...
    (S{4}') (S{5}') (S{6}') (S{7}')*theta(end-1)/4 (S{8}')*theta(end)/4])
legend({'1','2','3','4','5','6','7','8'})