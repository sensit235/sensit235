function [mnt_out, sdt_out] = sensit_morris(r,n,t,x0_min,x0_max,p_min,p_max,model,varargin)
%% run_morris.m
%
% This function computes the sensitivity measures for a given model using
% the Morris method. The sensitivity measures are the mean and standard
% deviation of the elementary effects associated with a parameter. The
% default behaviour is to normalise the sensitivity measures according to
% the standard definition. Optionally, the normalisation can be removed
% entirely or changed to that used for relative sensitivity functions,
% namely:
%
% $$\hat{\theta} = \frac{dy}{d\theta}\frac{\theta_0}{y}.$$
%
% Args:
%
% * |r| - # of random orientations (random starting points)
% * |n| - # number of divisions of parameter space
% * |t| - vector of times at which to compute sensitivity
% * |p_min| - vector of parameter min bound
% * |p_max| - vector of parameter max bound
% * |model| - a handle to the model RHS
% * |nls| - specify normalisation if required ('none','rsf')



% test for different normalisation
if isempty(varargin)
    nls = 'def';
else
    nls = varargin{1};
end

%%
% Add initial conditions to parameters to determine sensitivity and ensure
% the total number of initial conditions + parameters is even and generate
% experiments.

np = length(p_min);
ni = length(x0_min);
p_min = [p_min x0_min];
p_max = [p_max x0_max];

if ~(mod(length(p_min),2)==0) % test if not even
    p_min = [p_min 0.0];
    p_max = [p_max 1.0];
    p_even = false;
else
    p_even = true;
end

A = Generate_Experiment(length(p_min),n,r);

%%
% Loop through experiments aggregating results

y = {};
for i=1:length(A)
    
    % initialise results array
    r = zeros(size(A{i},1),length(t),ni);
    
    % loop through individual runs
    for j=1:size(A{i},1)
        
        % scale active parameters
        p_=scale_parameters(p_min,p_max,A{i}(j,:));
        
        % run model
        options = odeset('RelTol',1e-6,'AbsTol',1e-6);
        morris_model = @(t,x) model(t,x,p_(1:np));
        [t,x] = ode15s(morris_model,t,p_(np+1:np+ni),options);
        
        % save results
        
%         r(j,:) = x(:,1);
        r(j,:,:) = reshape(x(:,:),1,length(t),ni);
        
    end
    
    y{i} = r;
    
end

%%
% Prepare data structures for processing results

est={};
rst={};

% loop through experiments
for i=1:length(A)

    % loop through times
    for j=1:size(y{i},2)
        est(j,i) = {A{i}};
        
        % loop through state variables
        for k=1:ni
            rst(j,i,k) = {y{i}(:,j,k)'};
        end

    end

end

%%
% Compute standard Morris method mean and standard deviation

mnt = zeros(length(p_min),size(est,1),ni);
sdt = mnt;
for k=1:ni;
    for i=1:size(est,1)
        es=[est(i,:)];
        rs=[rst(i,:,k)];
        [mnt(:,i,k),sdt(:,i,k)] = Process_Results(es,rs);
    end
end

%%
% Remove normalisation and apply new if requested
est_orig = est;
if strcmp(nls,'none') | strcmp(nls,'rsf')
    for kp=1:ni
        ee = est;
        for j=1:size(est,2) % experiments
            for i=1:size(est,1) % times
                % remove normalisations
                for k=1:size(est{i,j},1)
    %                 pc=scale_parameters(theta,est{i,j}(k,:),a);
                    pc=scale_parameters(p_min,p_max,est{i,j}(k,:));
                    est{i,j}(k,:) = pc;
                end

                % apply new normalisation
                [ee{i,j},sd] = Process_Results(est(i,j),rst(i,j,kp));
                [l,~,~] = find(diff(est{i,j}));
                for k=1:size(ee{i,j},1)
                    ee{i,j}(k) = ee{i,j}(k).*est{i,j}(l(k),k)./rst{i,j,kp}(l(k));
                end
            end
        end

        % output with normalisation removed
        if strcmp(nls,'none')
            for i=1:size(est,1)
                es=est(i,:);
                rs=rst(i,:,kp);
                [mnt(:,i,kp),sdt(:,i,kp)] = Process_Results(es,rs);
            end
        end

        % output with relative sensitivity normalisation
        if strcmp(nls,'rsf')
            for i=1:size(est,1) % times
                mnt(:,i,kp) = mean(cell2mat(ee(i,:)),2);
                sdt(:,i,kp) = std(cell2mat(ee(i,:)),0,2);
            end
        end
        
        % reset est
        est = est_orig;
    end
end

%% Post processing
%
% Strip off extra parameter if odd number of parameters

if ~p_even
    mnt = mnt(1:end-1,:,:);
    sdt = sdt(1:end-1,:,:);
end

% construct array of sensitivities
mnt_out=[];
sdt_out=[];
for k=1:ni
    mnt_out = [mnt_out; mnt(:,:,k)];
    sdt_out = [sdt_out; sdt(:,:,k)];
end

% append mean solutions
solns = permute(mean(cat(1,y{:})),[3 2 1]);
mnt_out = [solns;mnt_out];

end