function [modelstrct,z,u,eta,models] = timing_var_BIC(X,Qvc,varargin)
% [modelstrct,z,u,eta,W,psi,sigma,iter,loglk] = timing_var_BIC(X,Qvc,varargin)
%
% MATLAB S3: Returns a structure (and parameter sets), containing model information associated
% with the lowest Bayesian Information Criteriorn (BIC), as described in Glaze &
% Troyer (2012).
%
% NOTE: In manuscript Sigma and Psi refer to component covariance
% matrices, while for effiency in the code, sigma and psi simply represent the diagonals
% of those matrices
%
%
% FUNCTION ARGUMENTS
% X: NxK matrix of data, N = # of samples, K is # of dimensions (i.e. unique time
% intervals per sequence)
%
% Qvc: Vector of dimensions in the global variable (BIC computed for each).
%
%
% OPTIONAL ARGUMENTS
%
% initN: Total # initial parameter estimates to run algorithm from in 1st
% stage in timing_var.m function (default=100).
%
% loglkN: # of initial parameters to choose associated with top loglkN
% log-likelihoods in 2nd stage in timing_var function (default=5).
%
% maxiter: Maximum number of iterations allowed for each EM run
% (default=10000).
%
% crit_loglk: Convergence criterion for log-likelihood in 1st stage in
% timing_var.m function (default=0.1%).
%
% crit_param: Convergence criterion for parameters in 2nd stage in
% timing_var.m function (default=0.01%).
%
% waitbaropt: For timing_var function, if set to 1, gives a waitbar 
% updating user on progress of code, if 0 no waitbar (default=1).
%
% psiopt: If psiopt=1, independent parameters are used (along with others) to
% evaluate convergence in 2nd stage. Default is 1.
%
% Optional arguments can be specified in any combination or order but
% require strings preceding each (as in other Matlab functions). For example:
%
% >> modelstrct = timing_var_BIC(X,[1 2 3],'initN',200,'loglkN',20);
%
% will compute model information for global dimension = 1 through 3 and
% return the model associated with the lowest BIC. Here, for each global dimension, 
% 200 initial conditions are used in the 1st stage, and the top 20 associated with the highest
% log-likelihoods chosen for the 2nd stage, with all other criteria set
% by default.
%
%
% FUNCTION RETURNS
% modelstrct: Structure containing model information associated with lowest
% BIC. Fields contain the following:
% 
% modelstrct.W: global weight matrix
%
% modelstrct.psi = independent variances
%
% modelstrct.sigma = jitter variances
%
% modelstrct.iter = total # iterations required for model to converge
%
% modelstrct.loglk = log-likelihood of the model
%
% modelstrct.Qvc = vector of global dimensions
%
% modelstrct.BIC = vector of BIC values computed for each dimension in Qvc
%
% modelstrct.SRMR = standardized root mean-squared residual for model
%
%
% z: NxQ matrix with global latent variables.
%
% u: NxK-1 matrix of jitter latent variables.
%
% eta: NxK matrix of independent latent variables.
%
% models: Array of models generated for all values in Qvc.
%
%
% AUTHOR: Chris Glaze (cmglaze@gmail.com)


% get user-secified arguments and assign values to the remaining
allargs = {'initN','loglkN','maxiter','crit_loglk','crit_param','waitbaropt','psiopt'};
args2timing_vars(allargs,varargin);

[N,K] = size(X);

modelnm = length(Qvc);
loglkvc = zeros(1,modelnm);

% Run timing_var.m for global dimension # specified in Qvc.
for modelind = 1:length(Qvc)  
    
    modelind
    
    models{modelind} = timing_var(X,Qvc(modelind),...
        'initN',initN,...
        'loglkN',loglkN,...
        'maxiter',maxiter,...
        'crit_loglk',crit_loglk,...
        'crit_param',crit_param,...
        'waitbaropt',waitbaropt,...
        'psiopt',psiopt);
    
    loglkvc(modelind) = models{modelind}.loglk;

end

% Compute # free parameters and BIC for each dimension #.
k = (Qvc + 2)*K - 1;
BIC = -2 * loglkvc + log(N) * k;

% Find minimum BIC and assign model information to modelstrct.
[minbic,modelind] = min(BIC);

modelstrct = models{modelind};
modelstrct.Qvc = Qvc;
modelstrct.BIC = BIC;

% Compute SRMR for model as described in main manuscript. 
R = corrcoef(X);
D = -diff(eye(K))';
Shat = modelstrct.W*modelstrct.W' + D*diag(modelstrct.sigma)*D' + diag(modelstrct.psi);
stdX = std(X);
Rhat = Shat ./ (stdX'*stdX);

% Indexing to unique entries (correlation and covariance matrices are
% symmetric, etc.).
u_inds = find(tril(ones(K)));

modelstrct.SRMR = sqrt(mean((R(u_inds)-Rhat(u_inds)).^2));

% Compute latent variables if required by user.
if nargout > 1
    [z,u,eta] = timing_latvar(X,modelstrct.W,modelstrct.sigma,modelstrct.psi);
end