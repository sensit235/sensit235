%% run_morris
%
% run the morris method to obtain sensitivity coefficients

function [] = run_morris()

% load up baseline parameters
parameters
p.dummy=0; % make up to even number of parameters

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

% function for scaling experiments to actual values based on %
function [p_]=scale_parameters(p,x,a)

    m=0; % keep track of how many parameters are changed
    p_=p;
    fn=fieldnames(p);
    for i=1:length(fn)
        if a(i)
            m=m+1;
            eval(['p_.' fn{i} '=p.' fn{i} '*(0.75 + 0.5*' num2str(x(m)) ');'])
        end
    end

end

% % function for scaling experiments to actual values based on ranges
% function [p_]=scale_parameters(p,pmin,pmax,x,a)
% 
%     m=0; % keep track of how many parameters are changed
%     p_=p;
%     fn=fieldnames(p);
%     for i=1:length(fn)
%         if a(i)
%             m=m+1;
%             eval(['p_.' fn{i} '=(pmax.' fn{i} '-pmin.' fn{i} ')*' num2str(x(m)) '+ pmin.' fn{i} ';'])
%         end
%     end
% 
% end