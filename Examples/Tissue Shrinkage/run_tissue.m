%% tissue.m
%
% This script computes the sensitivity of the tissue contraction model.
% This model has an observation function, which has been included by adding
% an extra equation to the system of ODEs.

close all
clear all
addpath ../../Methods/TSF/
addpath ../../Methods/Morris' Method'/


%
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
theta = [0.0100   -3.2798    0.0061   -2.3690...
    0.0040   -1.4710 0.00003964   0.00002701].*1e4;

% initial condition
x0 = [1 0 0 0]; % [N U D]

%% Run the tsf sensitivity analysis.
%
% The vector of parameters under consideration is constructed and the
% |sensit_tsf| function called to compute the sensitivities

T =@(t) (358-(273+37))*heaviside(-t+120)+273+37;
tsf_tissue_RHS =@(t,x,p) tissue_RHS(t,x,p,T);

tspan = [linspace(0,1,100) linspace(1.01,1500,100)]';
[t y] = sensit_tsf(tspan,tsf_tissue_RHS,theta,x0);

%% Compute sensitivities of observation function and plot.
ob1 = [-100,-100*theta(end-1),-100*theta(end)];

% parameters in system of ODEs
S = {};
for i=5:10
    S = [S ob1*(y(:,i:8:27)')];
end

% parameters in observation function
S = [S -100*(y(:,2)')];
S = [S -100*(y(:,3)')];

% compare direct computation of xi versus solution of 4th equation
figure
plot(tspan,100+ob1*([y(:,1) y(:,2) y(:,3)]'),tspan,y(:,4),'o')

% compare direct computation of sensitivity versus 4th equation
figure
plot(tspan,[(S{1}')*theta(1) (S{2}')*theta(2) (S{3}')*theta(3)...
    (S{4}')*theta(4) (S{5}')*theta(5) (S{6}')*theta(6)...
    (S{7}')*theta(7) (S{8}')*theta(8)])
hold on
plot(tspan,y(:,29:36)*diag(theta),'blackx')

%% Run the morris sensitivity analysis.
%
% In order to run the morris method the parameters must have an initial
% value > 0. At the moment the code with fail to run if this is not the
% case.

% amend x0
x0 = [0.999 0.0005 0.0005 100*(1-(0.999+theta(7)*0.0005+theta(8)*0.0005))];

pcg = 1e-1;
theta_min = (1 + pcg/100).*theta;
theta_max = (1 - pcg/100).*theta;
x0_min = (1 + pcg/100).*x0;
x0_max = (1 - pcg/100).*x0;

T =@(t) (358-(273+37))*heaviside(-t+120)+273+37;
morris_tissue_RHS =@(t,x,p) tissue_RHS(t,x,p,T);

tspan = [linspace(0,1,100) linspace(1.01,1500,100)]';

[mnt1 sdt1] = ...
    sensit_morris(4,4,tspan,x0_min,x0_max,theta_min,theta_max,morris_tissue_RHS,'rsf');

%% Compute sensitivities of observation function and plot.
ob1 = [-100,-100*theta(end-1),-100*theta(end)];

% parameters in system of ODEs
S = {};
for i=5:12
    S = [S ob1*mnt1(i:12:36,:)];
end

figure
plot(tspan,100+ob1*mnt1(1:3,:),tspan,y(:,4),'o')

figure
plot(tspan,[(S{1}')*theta(1) (S{2}')*theta(2) (S{3}')*theta(3)...
    (S{4}')*theta(4) (S{5}')*theta(5) (S{6}')*theta(6)])
legend({'1','2','3','4','5','6'})
hold on
plot(tspan,(mnt1(41:48,:)')*diag(theta),'blackx')

%% display useful results

% tsf versus morris
figure
plot(tspan,mnt1(41:48,:)',...
    tspan,y(:,29:36)./kron(ones(1,8),y(:,4))*diag(theta),'blackx')
legend({'1','2','3','4','5','6','7','8'})

% solution / state variables / sensitivities
figure
subplot(3,1,1)
plot(tspan,y(:,4))
subplot(3,1,2)
plot(tspan,[y(:,1) y(:,2) y(:,3)])
subplot(3,1,3)
plot(tspan,mnt1(41:48,:)')
legend({'1','2','3','4','5','6','7','8'})
