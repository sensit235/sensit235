%% run_gen_exp
%
% generate the experiments and save to file for run_morris
%
% args:
%    np - # parameters
%    nx - # experiments
%    nd - # of divisions per parameter

function [] = run_gen_exp(np, nx, nd)

    B = Generate_Experiment(np,nd,nx);

    for i=1:length(B)
        A = B{i};
        eval(['save experiment-' num2str(i) '.mat A'])
    end

end