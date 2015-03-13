%% benchmark.m
%
% a set of benchmarks for the morris method

close all
clear all

dx=0.01; % for plotting
n = 10; % number of random orientations (higher is better accuracy)

%% run test function 1
%
% a linear function with 4 parameters
%    y = a*x1 + b*x2 + c*x3 + d*x4

A = Generate_Experiment(4,4,n);

for i=1:n
    r{i} = test_function_1(A{i}');
end

[mn sd] = Process_Results(A,r)

subplot(3,1,1)
plot(mn,sd,'s','LineWidth',2)
hold on
plot([90.050 71.045 41.470 20.825],[31.08 28.86 17.75 11.54],'or')
xlabel('Mean')
ylabel('Standard Deviation')
text(abs(mn)+dx.*norm(mn),sd+dx*norm(sd),cellstr(num2str([1:length(mn)]')))
legend('morris method','exact','Location','SouthEast')

%% run test function 2 with same experimental design
%
% second order interations included

for i=1:n
    r{i} = test_function_2(A{i}');
end

[mn sd] = Process_Results(A,r);

subplot(3,1,2)
plot(mn,sd,'s','LineWidth',2)
hold on
plot([0.05 0.59 10 0.21],[0 0 0 0],'or')
xlabel('Mean')
ylabel('Standard Deviation')
text(abs(mn)+dx.*norm(mn),sd+dx*norm(sd),cellstr(num2str([1:length(mn)]')))
legend('morris method','exact','Location','SouthEast')

%% run test function 3
%
% sobol g-function
A = Generate_Experiment(8,4,n);

for i=1:n
    r{i} = test_function_3(A{i}');
end

[mn sd] = Process_Results(A,r);

subplot(3,1,3)
plot(mn,sd,'s','LineWidth',2)
xlabel('Mean')
ylabel('Standard Deviation')
text(mn+dx.*norm(mn),sd+dx*norm(sd),cellstr(num2str([1:length(mn)]')))
legend('morris method','exact','Location','SouthEast')