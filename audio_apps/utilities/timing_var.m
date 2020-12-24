function [modelstrct,z,u,eta] = timing_var(X,Q,varargin)
% [modelstrct,z,u,eta] = timing_var(X,Q,varargin)
% 
% MATLAB S2: Returns parameters, latent variable estimates from the timing variability model, 
% running the EM algorithm from multiple initial conditions and picking the parameter set associated
% with the highest log-likelihood.
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
% Q: Dimensionality of global variable.
%
%
% OPTIONAL ARGUMENTS
%
% initN: Total # initial parameter estimates to run algorithm from in 1st
% stage (default=100).
%
% loglkN: # of initial parameters to choose associated with top loglkN
% log-likelihoods in 2nd stage (default=5).
%
% maxiter: Maximum number of iterations allowed for each EM run
% (default=10000).
%
% crit_loglk: Convergence criterion for log-likelihood in 1st stage
% (default=0.1%).
%
% crit_param: Convergence criterion for parameters in 2nd stage
% (default=0.01%).
%
% waitbaropt: If set to 1, gives a waitbar updating user on progress of
% code, if 0 no waitbar (default=1).
%
% psiopt: If psiopt=1, independent parameters are used (along with others) to
% evaluate convergence in 2nd stage. Default is 1.
%
% Optional arguments can be specified in any combination or order but
% require strings preceding each (as in other Matlab functions). For example:
%
% >> [modelstrct,z,u,eta] = timing_var(X,2,'initN',200,'loglkN',20,'crit_param',0.1);
%
% will return estimates for data X, global dimension=2, with 200 initial
% conditions in the 1st stage, the top 20 associated with the highest
% log-likelihoods chosen for the 2nd stage, a convergence criterion
% in the 2nd stage of 0.1%, and all other arguments set by default.
%
%
% FUNCTION RETURNS
%
% modelstrct: Structure with model parameter estimates and information.
% Fields contain the following:
%
% modelstrct.W: KxQ matrix of global weights, K is # of dimensions in data set, Q is #
% of global dimensions.
%
% modelstrct.sigma: (K-1)x1 vector of jitter variances.
%
% modelstrct.psi: Kx1 vector of independent variances.

% modelstrct.loglk: log likelihood of the model.
%
% modelstrct.iter: total # of iterations needed for convergence for parameter set
% chosen.
%
%
% z: NxQ matrix with global latent variables.
%
% u: NxK-1 matrix of jitter latent variables.
%
% eta: NxK matrix of independent latent variables.
%
%
% AUTHOR: Chris Glaze (cmglaze@gmail.com)


% get user-secified arguments and assign values to the remaining
allargs = {'initN','loglkN','maxiter','crit_loglk','crit_param','waitbaropt','psiopt'};
args2timing_vars(allargs,varargin);

[N,K] = size(X);

% Initialize matrices and vectors to save parameters, etc. associated with 
% top loglkN log-likelihoods (a procedure is used to avoid having to store
% information for all initN models).
Wmat = zeros(K,Q,loglkN);
sigmamat = zeros(K-1,loglkN);
psimat = zeros(K,loglkN);
loglkvc = zeros(1,loglkN);
itervc = loglkvc;

init_ind = 0;
n_fail = 0;

if waitbaropt
    h = waitbar(0/initN,'Calculating FA parameters');
end


% STAGE 1: Loop criteria require EM to run until initN parameter sets have converged 
% or EM has failed to converge 500 times (in which case code is stopped).

while init_ind < initN && n_fail < 500
    
    % Run EM with log-likelihood convergence. Note that W is not rotated in
    % this stage so final W estimate can be used as initial condition in
    % 2nd stage.
    modelstrct = timing_var_EM(X,Q,...
        'loglkflg',1,...
        'conv_crit',crit_loglk,...
        'rotopt',0,...
        'psiopt',psiopt);
    
    % If convergence is achieved (i.e. iterations < maximum allowed), and
    % no degeneracy is detected, save model information.
    if modelstrct.iter < maxiter && modelstrct.degflg
        init_ind = init_ind + 1;
        
        % If fewer than loglkN parameter sets have converged, continue to
        % populate matrices and vectors that will provide information for
        % 2nd stage.
        if init_ind <= loglkN
            Wmat(:,:,init_ind) = modelstrct.W;
            sigmamat(:,init_ind) = modelstrct.sigma;
            psimat(:,init_ind) = modelstrct.psi;
            loglkvc(init_ind) = modelstrct.loglk;
            itervc(init_ind) = modelstrct.iter;
            
            
        % Otherwise, if the current log-likelihood is higher than the
        % lowest already stored, replace that model information with the current.
        else
            [loglk_min,minind] = min(loglkvc);
            if modelstrct.loglk > loglk_min
                Wmat(:,:,minind) = modelstrct.W;
                sigmamat(:,minind) = modelstrct.sigma;
                psimat(:,minind) = modelstrct.psi;
                loglkvc(minind) =  modelstrct.loglk;
                itervc(minind) = modelstrct.iter;
                
            end
            
        end
        
        if waitbaropt
            h = waitbar(init_ind/initN,h,'Calculating FA parameters');
        end
        
    else
        n_fail = n_fail + 1;
    end
    
    % If > 20 failures to converge have occured, start warning the user.
    if n_fail > 20
        warning(['Warning: ' num2str(n_fail) ' failures to converge.'])
    end
    
    clear modelstrct
    
end

if waitbaropt
    close(h)
end


% STAGE 2 (only if EM has run successfully initN times): Now take information 
% from top loglkN models and run EM until parameter convergence is achieved, 
% use parameters associated with highest log-likelihood resulting from that convergence.

if init_ind == initN
    loglkprev = -inf;
    
    if waitbaropt
        h = waitbar(0/loglkN,['Calculating top ' num2str(loglkN) ' model parameters']);
    end
    
    for loglk_ind = 1:loglkN
        
        % Run EM using parameter convergence. Final parameter
        % estimates from 1st stage are used as initial conditions so EM can
        % simply begin where it left off.
        modelstrct_tmp = timing_var_EM(X,Q,...
            'loglkflg',0,...
            'conv_crit',crit_param,...
            'rotopt',1,...
            'psiopt',psiopt,...
            'W0',squeeze(Wmat(:,:,loglk_ind)),...
            'psi0',psimat(:,loglk_ind),...
            'sigma0',sigmamat(:,loglk_ind));
        
        % If log-likelihood is higher than last, convergence was achieved 
        % and no degeneracy is detected, assign final variables to that model
        if modelstrct_tmp.loglk > loglkprev && modelstrct_tmp.iter < maxiter && modelstrct_tmp.degflg
            modelstrct = modelstrct_tmp;
            modelstrct.iter = modelstrct.iter + itervc(loglk_ind);
            
            loglkprev = modelstrct.loglk;
        end
        
        if waitbaropt
            h = waitbar(loglk_ind/loglkN,h,['Calculating top ' num2str(loglkN) ' model parameters']);
        end
        
    end
    
    if waitbaropt
        close(h)
    end
    
end

% If W does not exist it means that no parameter sets converged in the 2nd
% stage.
if ~exist('modelstrct')
    error('Error: no estimates converged')
end

modelstrct = rmfield(modelstrct,'degflg');

% Compute latent variables if required by output.
if nargout > 1
    [z,u,eta] = timing_latvar(X,modelstrct.W,modelstrct.sigma,modelstrct.psi);
end
