%% run_methanogenesis
%
% This script serves as an example of how to run methanogenesis.

% generate parameters
parameters

% data
data_1=[0 2 4 5 6 10 12 14 16 18];
data_2=[3.68 2.54 1.90 1.85 1.61 1.31 1.28 1.17 1.02 0.97]*1e-3;


% set initial values and time
tspan = [0 20*24*60*60];
x0 = [3.5e-3 45.2e-3 1e-6 0.009]'; % [molal molal molal g/kg]

% obtain solution
[t x] = methanogenesis(tspan,p,x0);

% plot
figure(1);
subplot(411)
    plot(t/60/60/24,x(:,1))
    title('m_{Ac}')
    xlabel('time (days)')
    ylabel('conc (molal)')
    hold on
    scatter(data_1,data_2)
subplot(412)
    plot(t/60/60/24,x(:,2))
    title('m_{HCO_3^-}')
    xlabel('time (days)')
    ylabel('conc (molal)')
subplot(413)
    plot(t/60/60/24,x(:,3))
    title('m_{CH_4}')
    xlabel('time (days)')
    ylabel('conc (molal)')
subplot(414)
    plot(t/60/60/24,x(:,4))
    title('[X]')
    xlabel('time (days)')
    ylabel('conc (molal)')