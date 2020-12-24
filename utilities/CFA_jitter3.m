function [W,sigma,psi,phi,Z,iter,degflg,z,u,q,j] = CFA_jitter3(X,Q,m,mz,W0,sigma0,psi0,phi0,Z0)
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


% maximum number of iterations
maxiter = 1000;

[n,d] = size(X);

X = X - repmat(mean(X),n,1);
% X = X ./ repmat(std(X),n,1);

% calculate jitter weight matrix in D, where multiplication by D is
% equivalent to differencing
D = inv(tril(ones(d)));
A = D;
Ainv = inv(A);

Qd = size(Q,2);

% calculate biased covariance matrix in C, variances in vr
C = X'*X/n;
vr = diag(C);

if nargin < 5
    
    % if no initial conditions are given, choose them at random by scaling
    % standard deviation or variances by a random number in [0,1]

    W0 = repmat(sqrt(vr),1,m).*rand(d,m);
    Z0 = repmat(sqrt(vr),1,mz).*rand(d,mz);
    

    sigma0 = diag(vr.*rand(d,1));
    psi0 = diag(vr.*rand(d,1));

    [rows,cols] = find(Q);
    rowsu = unique(rows);

    phi0 = diag(vr(rowsu(1:Qd)).*rand(Qd,1));

    %     sigma0 = diag(diag(C));
    %     psi0 = sigma0;
    %     phi0 = diag(rand(1,Qd));


end

% initalized parameter estimates

W = W0;
Z = Z0;
sigma = sigma0;
psi = psi0;
phi = phi0;

if nargout > 6
    logp = zeros(maxiter,1);
end

iter = 1;

sigmamat = zeros(maxiter,d);
psimat = sigmamat;
wmat = zeros(maxiter,d,m);
phimat = zeros(maxiter,max(Qd,1));
zmat = zeros(maxiter,d,mz);

sigmamat(iter,:) = diag(sigma);
psimat(iter,:) = diag(psi);
wmat(iter,:,:) = W;
zmat(iter,:,:) = Z;


% throughout Chat is the estimated covariance matrix, log-likelihood
% follows Bishop (2006)
Chat = W*W' + D*sigma*D' + Q*phi*Q' + psi + A*Z*Z'*A';
logp(iter) = -n/2*(d*log(2*pi) + log(det(Chat)) + trace(inv(Chat)*C));

if Qd > 0
    phimat(iter,:) = diag(phi);
end

detsigma = det(sigma);
detpsi = det(psi);
detphi = det(phi);

% deltaflg = 1 indicates lack of convergence
deltaflg = 1;

% degflg indicates degeneracy of at least 1 non-global matrix (expected global
% variance = 1 always)
degflg = detsigma > 0 & detpsi > 0 & detphi > 0;


while iter < maxiter && degflg && deltaflg

    % EM algorithm for estimating parameters as described in Glaze (2008),
    % based in EM for traditional factor analysis in Bishop (2006)
    
    invpsi = inv(psi);
   
    WDQZ = [W D Q A*Z];

    C_zuqj = [eye(m) zeros(m,d+Qd+mz);zeros(d,m) sigma zeros(d,Qd+mz);zeros(Qd,d+m) phi zeros(Qd,mz); zeros(mz,d+m+Qd) eye(mz)];   
    G_zuqj = inv(inv(C_zuqj) + WDQZ'*invpsi*WDQZ);
    
    zuqj = X*invpsi*WDQZ*G_zuqj;
    zuqj2sm = zuqj'*zuqj + n*G_zuqj;
    
    sigma = diag(diag(zuqj2sm(m+1:m+d,m+1:m+d))) / n;
    phi = diag(diag(zuqj2sm(m+d+1:m+d+Qd,m+d+1:m+d+Qd))) / n;

    WDQZ = (X'*zuqj * inv(zuqj2sm));
    
    W = WDQZ(:,1:m);
    Z = Ainv*WDQZ(:,end-mz+1:end);

    psifull = C + WDQZ*zuqj2sm*WDQZ'/n - 2*X'*zuqj*WDQZ'/n;
    psi = diag(diag((psifull)));
    
    
    
    detsigma = det(sigma);
    detpsi = det(psi);
    detphi = det(phi);

    iter = iter + 1;

    sigmamat(iter,:) = diag(sigma);
    psimat(iter,:) = diag(psi);
    wmat(iter,:,:) = W;
    
    if nargout > 9
        Chat = W*W' + D*sigma*D' + Q*phi*Q' + psi + A*Z*Z'*A';
        logp(iter) = -n/2*(d*log(2*pi) + log(det(Chat)) + trace(inv(Chat)*C));
    end

    if Qd > 0
        phimat(iter,:) = diag(phi);
    end

    
    % settle on estimates if all variance parameters have changed by < 0.01
    % squared-units from the pervious iteration and global weight estimates
    % have changed by < 0.1 units
    
    if max(abs(diff(sigmamat(iter-1:iter,:)))) < 0.01 && max(abs(diff(psimat(iter-1:iter,:)))) < 0.01 ...
            &&  max(max(abs(diff(wmat(iter-1:iter,:,:))))) < 0.1 &&  max(abs(diff(phimat(iter-1:iter,:)))) < 0.01 &&  max(max(abs(diff(zmat(iter-1:iter,:,:))))) < 0.1
        deltaflg = 0;
    end

    degflg = detsigma > 0 & detpsi > 0 & detphi > 0;

end

psimat = psimat(1:iter,:);
sigmamat = sigmamat(1:iter,:);
wmat = wmat(1:iter,:,:);
phimat = phimat(1:iter,:);
zmat = zmat(1:iter,:,:);

if nargout > 9
    logp = logp(1:iter);
end

z = zuqj(:,1:m);
u = zuqj(:,m+1:m+d);
q = zuqj(:,m+d+1:m+d+Qd);
j = zuqj(:,m+d+Qd+1);
