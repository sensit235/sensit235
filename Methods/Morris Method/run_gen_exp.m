%% run_gen_exp
%
% generate the experiments and save to file for run_morris

np = 6; % # parameters
nx = 10; % # experiments
nd = 4; % number of divisions of parameter

B = Generate_Experiment(np,nd,nx);

for i=1:length(B)
    A = B{i};
    eval(['save experiment-' num2str(i) '.mat A'])
end