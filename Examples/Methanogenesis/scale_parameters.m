%% scale_parameters.m
%
%

% function for scaling experiments to actual values based on %
function [p_]=scale_parameters(p,x,a)

    pcg = 0.0001; % chosen fraction to vary params by
    m=0; % keep track of how many parameters are changed
    p_=p;
    fn=fieldnames(p);
    for i=1:length(fn)
        if a(i)
            m=m+1;
%             eval(['p_.' fn{i} '=p.' fn{i} '*(0.75 + 0.5*' num2str(x(m)) ');'])
            eval(['p_.' fn{i} '=p.' fn{i} '*(1.0 - ' num2str(pcg) ' + 2*' num2str(pcg) '*' num2str(x(m)) ');'])
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