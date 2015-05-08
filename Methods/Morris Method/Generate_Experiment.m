%% Generate_Experiment.m
%
% This function will generate a set of computational experiments that are
% used in the Morris method to compute the sensitivity measures. Each
% experiment is a $(p+1)\times p$ matrix with each row representing a
% single sample of the parameter space.
%
% Args:
%
% * _k_ = # parameters
% * _p_ = # divisions for each parameter
% * _r_ = # random orientations
%
% Returns:
%
% * _experiments_ = cell array of experiments

function experiments = Generate_Experiment(k,p,r)

if rem(k,2)
    error('Only valid for EVEN number of parameters')
end

% initialise design matrices
B = tril(ones(k+1,k),-1);
D_ = diag(randi(2,1,k)*2-3); % random diagonal
J = ones(k+1,k);

% generate random orientation
for i=1:r
    x_ = (randi(p-1,1,k)-1)./(p-1); % pick random parameter vector
    P_ = eye(k);
    P_ = P_(randperm(k),:); % random permutation matrix

    % random orientation of B
    experiments{i} = (J(:,1)*x_ + (1/(2*(p-1))).*((2*B-J)*D_+J))*P_;
end

end