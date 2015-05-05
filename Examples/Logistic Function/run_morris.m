%% run_morris
%
% run the morris method to obtain sensitivity coefficients

% args:
%   t - equivalent of tspan in ode15s - time range of solution
%   p - vector of parameters
%   a - boolean indicating if parameter is included

function [] = run_morris(t,p,a)

tspan = t;

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
        
        % run model
        [t x] = logistic(tspan,p_,p_(3));
        
        % save results
        r(j,:) = x(:,1);
        
    end
    
    eval(['save results-' file_list(i).name ' r'])
    
end

end