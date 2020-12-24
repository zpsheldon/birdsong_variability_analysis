function [sigma,psi,sigma2,iter,degflg,u,j,logp,sigmamat] = CFA_jitteronly2(X,J,sigma0,psi0,sigma20)

% maximum number of iterations
maxiter = 1000;

[n,d] = size(X);

X = X - repmat(mean(X),n,1);
% X = X ./ repmat(std(X),n,1);

% calculate jitter weight matrix in D, where multiplication by D is
% equivalent to differencing
D = inv(tril(ones(d)));

% calculate biased covariance matrix in C, variances in vr
C = X'*X/n;
vr = diag(C);

if nargin < 3
    
    % if no initial conditions are given, choose them at random by scaling
    % standard deviation or variances by a random number in [0,1]

    sigma0 = diag(vr.*rand(d,1));
    sigma20 = mean(vr) * rand;
    psi0 = diag(vr.*rand(d,1));

    %     sigma0 = diag(diag(C));
    %     psi0 = sigma0;
    %     phi0 = diag(rand(1,Qd));


end

% initalized parameter estimates

sigma = sigma0;
sigma2 = sigma20;
psi = psi0;

if nargout > 6
    logp = zeros(maxiter,1);
end

iter = 1;

sigmamat = zeros(maxiter,d);
psimat = sigmamat;
sigma2mat = zeros(maxiter,1);

sigmamat(iter,:) = diag(sigma);
psimat(iter,:) = diag(psi);
sigma2mat(iter) = sigma2;

% throughout Chat is the estimated covariance matrix, log-likelihood
% follows Bishop (2006)
Chat = D*sigma*D' + psi + J*sigma2*J';
logp(iter) = -n/2*(d*log(2*pi) + log(det(Chat)) + trace(inv(Chat)*C));

detsigma = det(sigma);
detpsi = det(psi);

% deltaflg = 1 indicates lack of convergence
deltaflg = 1;

% degflg indicates degeneracy of at least 1 non-global matrix (expected global
% variance = 1 always)
degflg = detsigma > 0 & detpsi > 0 & sigma2 > eps;

DJ = [D J];

while iter < maxiter && degflg && deltaflg

    % EM algorithm for estimating parameters as described in Glaze (2008),
    % based in EM for traditional factor analysis in Bishop (2006)
    
    invpsi = inv(psi);

    C_uj = [sigma zeros(d,1);zeros(1,d) sigma2];
    G_uj = inv(inv(C_uj) + DJ'*invpsi*DJ);

    uj = X*invpsi*DJ*G_uj;
    uj2sm = uj'*uj + n*G_uj;

    sigma = diag(diag(uj2sm(1:d,1:d))) / n;
    sigma2 = uj2sm(d+1,d+1) / n;

    psifull = C + DJ*uj2sm*DJ'/n - 2*X'*uj*DJ'/n;
    psi = diag(diag((psifull)));


    detsigma = det(sigma);

    iter = iter + 1;

    sigmamat(iter,:) = diag(sigma);
    psimat(iter,:) = diag(psi);
    sigma2mat(iter) = sigma2;
    
    if nargout > 9
        Chat = D*sigma*D' + psi + J*sigma2*J';
        logp(iter) = -n/2*(d*log(2*pi) + log(det(Chat)) + trace(inv(Chat)*C));
    end

    
    % settle on estimates if all variance parameters have changed by < 0.01
    % squared-units from the pervious iteration and global weight estimates
    % have changed by < 0.1 units
    
    if max(abs(diff(sigmamat(iter-1:iter,:)))) < 0.01 && max(abs(diff(psimat(iter-1:iter,:)))) < 0.01  ...
            && max(abs(diff(sigma2mat(iter-1:iter,:)))) < 0.01
        deltaflg = 0;
    end

    degflg = detsigma > 0 & detpsi > 0 & sigma2 > eps;

end

psimat = psimat(1:iter,:);
sigmamat = sigmamat(1:iter,:);
sigma2mat = sigma2mat(1:iter);

if nargout > 9
    logp = logp(1:iter);
end

u = uj(:,1:d);
j = uj(:,d+1);

