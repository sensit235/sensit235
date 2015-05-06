%% run_morris
%
% run the morris method to obtain sensitivity coefficients

function [] = run_morris()

% load up baseline parameters
parameters

% set time
tspan = linspace(0,20*24*60*60,100);

% initial condition
x0 = [3.5e-3 45.2e-3 1e-6 0.009]'; % [molal molal molal g/kg]

% define active parameters (those being varied)
a = [false false true true false true true false true false true false];

% find files holding experiments
file_list=dir('experiment-*');
nf=length(file_list);

% loop through experiments
for i=1:nf
    
    eval(['load ' file_list(i).name]);
    
    % initialise results array
    r = zeros(size(A,1),length(tspan));
    
    % loop through individual runs
    for j=1:size(A,1)
        
        % scale active parameters
        p_=scale_parameters(p,A(j,:),a);
%         p_=scale_parameters(p,pmin,pmax,A(j,:),a);
        
        % run model
        [t x] = methanogenesis(tspan,p_,x0);
        
        % save results
        r(j,:) = x(:,1);
        
    end
    
    eval(['save results-' file_list(i).name ' r'])
    
end

end