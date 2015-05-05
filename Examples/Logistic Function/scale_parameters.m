%% scale_parameters.m
%
%

% function for scaling experiments to actual values based on %
function [p_]=scale_parameters(p,x,a)

    pcg = 0.0000001; % chosen fraction to vary params by
    m=0; % keep track of how many parameters are changed
    p_=p;
    for i=1:length(p)
        if a(i)
            m=m+1;
            p_(m) = p(m)*(1.0 - pcg + 2*pcg*x(m));
        end
    end

end