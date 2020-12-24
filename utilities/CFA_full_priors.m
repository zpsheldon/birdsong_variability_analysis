function [W,sigma,psi,phi,omega,iter,degflg,z,u,q,j,G_zuqj,logpvc,wmat,sigmamat,psimat,phimat,omegamat] = CFA_full_priors(X,Q,m,W0,sigma0,psi0,phi0,omega0,D,R,Wexp,conv_crit,logpflg)
% [W,sigma,psi,phi,iter,degflg,z,u,q,logp,wmat,sigmamat,psimat,phimat] = CFA_jitter(X,Q,m,W0,sigma0,psi0,phi0)
% returns parameters of the factor analysis model described in Glaze (2008)
%
%
% FUNCTION RETURNS
% W: dXm matrix of factor weights, d is # of dimensions in data set, m is #
% of factor weights.
%
% sigma: dX1 vector of jitter variances
%
% psi: dX1 vector of residual (independent) variances
%
% phi: pX1 matrix of repeated variances, p is # of unique dimensions that
% are repeated
%
% iter: number of iterations needed for convergence
% degflg: flag indicating whether routine was terminated due to a
% degenerate matrix
%
% z: nXm matrix with global factors
% u: nXd matrix of jitter factors
% q: nXp matrix of repeated independent factors
%
% logp: iterX1 vector the model log probability after each iteration
% wmat, sigmamat, psimat, phimat: iterXdXm, iterXd, iterXd, and iterXp
% matrices of parameter estimates after each iteration
%
% FUNCTION ARGUMENTS
% X: nXd matrix of data, n = # of samples, d is # of dimensions
%
% Q: a weight matrix for repeated noise, each row corresponds to a
% dimension of X, each column is a "repeat-factor". Q(i,j) = 1 if dimension
% i is influenced by factor j, 0 otherwise
%
% m is the number of global parameters to estimate
%
% W0, sigma0, psi0, phi0: optional arguments giving initial conditions for
% parameter estimates


if nargin < 13
    logpflg = 0;
end

% maximum number of iterations
maxiter = 100000;

[n,d] = size(X);

X = X - repmat(mean(X),n,1);
% X = X ./ repmat(std(X),n,1);

% calculate jitter weight matrix in D, where multiplication by D is
% equivalent to differencing

if nargin < 9
    D = inv(tril(ones(d)));
end

if nargin < 10
    R = D*Q;
end

if nargin < 12
    conv_crit = 0.01;
end


Dd = size(D,2);
Qd = size(Q,2);
Rd = size(R,2);

% calculate biased covariance matrix in C, variances in vr
C = X'*X/(n);
vr = diag(C);

if nargin < 4
    
    % if no initial conditions are given, choose them at random by scaling
    % standard deviation or variances by a random number in [0,1]
    
    W0 = repmat(sqrt(vr),1,m).*rand(d,m);
    
    psi0 = vr.*rand(d,1);
    
    sigma0 = diag(D'*diag(vr.*rand(d,1))*D ./ diag(diag(D'*D)));
    phi0 = diag(Q'*diag(vr.*rand(d,1))*Q ./ diag(diag(Q'*Q)));
    omega0 = diag(R'*diag(vr.*rand(d,1))*R ./ diag(diag(R'*R)));
    
end

% initalized parameter estimates

W = W0;
sigma = sigma0;
psi = psi0;
phi = phi0;
omega = omega0;

if nargout > 12 || logpflg
    logpvc = zeros(maxiter,1);
end

iter = 1;

if nargout > 13
    
    sigmamat = zeros(maxiter,max(Dd,1));
    psimat = zeros(maxiter,d);
    wmat = zeros(maxiter,d,max(m,1));
    phimat = zeros(maxiter,max(Qd,1));
    omegamat = zeros(maxiter,max(Rd,1));
    
    psimat(iter,:) = psi;
    
    if m > 0
        wmat(iter,:,:) = W;
    end
    
    if Qd > 0
        phimat(iter,:) = phi;
    end
    
    if Rd > 0
        omegamat(iter,:) = omega;
    end
    
    if Dd > 0
        sigmamat(iter,:) = sigma;
    end
    
end

% throughout Chat is the estimated covariance matrix, log-likelihood
% follows Bishop (2006)
Chat = W*W' + D*diag(sigma)*D' + Q*diag(phi)*Q' + diag(psi) + R*diag(omega)*R';
logpvc(iter) = -n/2*(d*log(2*pi) + log(det(Chat)) + trace(inv(Chat)*C)) -.5*sum((W(:)-Wexp(:)).^2)-sum(2*log(psi)+psi.^-2)-sum(2*log(sigma)+sigma.^-2)-sum(2*log(phi)+phi.^-2)-sum(2*log(omega)+omega.^-2);



detsigma = prod(sigma);
detpsi = prod(psi);
detphi = prod(phi);
detomega = prod(omega);

% deltaflg = 1 indicates lack of convergence
deltaflg = 1;

% degflg indicates degeneracy of at least 1 non-global matrix (expected global
% variance = 1 always)
degflg = detsigma > 0 & detpsi > 0 & detphi > 0 & detomega > 0;

DQR = [D Q R];
WDQR = [W DQR];

mones = ones(m,1);
C_zuqj = [mones;sigma;phi;omega];

I = eye(d);

while iter < maxiter && degflg && deltaflg
    
    % EM algorithm for estimating parameters as described in Glaze (2008),
    % based in EM for traditional factor analysis in Bishop (2006)
    
    invpsi = diag(1./psi);
    G_zuqj = inv(diag(1./C_zuqj) + WDQR'*invpsi*WDQR);
        
    % this step involves a lot of multiplication by 0 and can probably be made much
    % more efficient
    zuqj = X*invpsi*WDQR*G_zuqj;    
    
    zuqj2sm = zuqj'*zuqj + n*G_zuqj;
    
    sigma = diag((zuqj2sm(m+1:m+Dd,m+1:m+Dd))+1) / (n+1);
    phi = diag(zuqj2sm(m+Dd+1:m+Dd+Qd,m+Dd+1:m+Dd+Qd)+1) / (n+1);
    omega = diag(zuqj2sm(m+Dd+Qd+1:m+Dd+Qd+Rd,m+Dd+Qd+1:m+Dd+Qd+Rd)+1) / (n+1);
    
    %     invzuqj2sm = inv(zuqj2sm);
    %     W = X'*zuqj(:,1:m)*invzuqj2sm(1:m,1:m);
    %
    
    Xzuqj = X'*zuqj;
%     Xz = X'*zuqj(:,1:m);
%     invzuqj2sm = inv(zuqj2sm);
  
    invz = inv(zuqj2sm(1:m,1:m)+1);

    W_prev = W;
%     W = Xzuqj(:,1:m) * invzuqj2sm(1:m,1:m);
  
    W = (Xzuqj(:,1:m)-DQR*zuqj2sm(m+1:end,1:m)+Wexp') * invz;

    C_zuqj_prev = C_zuqj;
    C_zuqj = [mones;sigma;phi;omega];
    
    WDQR = [W DQR];
    
    psi_prev = psi;
    psi = diag(n*(C + WDQR*zuqj2sm*WDQR'/(n) - 2*Xzuqj*WDQR'/(n))+1)/(n+1);
    
    detsigma = prod(sigma);
    detpsi = prod(psi);
    detphi = prod(phi);
    detomega = prod(omega);
    
    iter = iter + 1;
    
    if nargout > 12 || logpflg
        Chat = W*W' + D*diag(sigma)*D' + Q*diag(phi)*Q' + diag(psi) + R*diag(omega)*R';
        logpvc(iter) = -n/2*(d*log(2*pi) + log(det(Chat)) + trace(inv(Chat)*C))-.5*sum((W(:)-Wexp(:)).^2)-sum(2*log(psi)+psi.^-2)-sum(2*log(sigma)+sigma.^-2)-sum(2*log(phi)+phi.^-2)-sum(2*log(omega)+omega.^-2);
    end
    
    if nargout > 13
    
        if Qd > 0
            phimat(iter,:) = phi;
        end
        
        if m > 0
            wmat(iter,:,:) = W;
        end
        
        if Dd > 0
            sigmamat(iter,:) = sigma;
        end
        
        if Rd > 0
            omegamat(iter,:) = omega;
        end
        
        psimat(iter,:) = psi;
    
    end
    
    % settle on estimates if all variance parameters have changed by <
    % conv_crit squared-units from the pervious iteration and global weight estimates
    % have changed by < sqrt(conv_crit) units
    
    if m == 0
        Wdiff = 0;
    else
        Wdiff = 100 * max(abs(1 - W(:) ./ W_prev(:)));
    end
    
    if isempty(C_zuqj)
        Cdiff = 0;
    else
        Cdiff = 100 * max(abs(1 - sqrt(C_zuqj) ./ sqrt(C_zuqj_prev)));
    end
    
    psidiff = 100 * max(abs(1 - sqrt(psi) ./ sqrt(psi_prev)));
    
    if logpflg
        deltaflg = ~(logpvc(iter) - logpvc(iter-1) > 0 && 100*abs((logpvc(iter) - logpvc(iter-1))/logpvc(iter-1)) < conv_crit);
    else
%         deltaflg = ~(Cdiff < conv_crit && Wdiff < sqrt(conv_crit));

        deltaflg = ~(psidiff < conv_crit && Cdiff < conv_crit && Wdiff < conv_crit);
    end
    degflg = detsigma > 0 & detpsi > 0 & detphi > 0 & detomega > 0;
    
end

if nargout > 13
    psimat = psimat(1:iter,:);
    sigmamat = sigmamat(1:iter,:);
    wmat = wmat(1:iter,:,:);
    phimat = phimat(1:iter,:);
    omegamat = omegamat(1:iter,:);    
end

if nargout > 13 || logpflg
    logpvc = logpvc(1:iter);
end

z = zuqj(:,1:m);
u = zuqj(:,m+1:m+Dd);
q = zuqj(:,m+Dd+1:m+Dd+Qd);
j = zuqj(:,m+Dd+Qd+1:m+Dd+Qd+Rd);

invpsi = diag(1./psi);
G_zuqj = inv(diag(1./C_zuqj) + WDQR'*invpsi*WDQR);
