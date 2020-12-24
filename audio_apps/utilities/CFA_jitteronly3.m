function [Z,psi,sigma,iter,degflg,u,j,logp,zmat,psimat,sigmamat] = CFA_jitteronly3(X,Z0,psi0,sigma0,mz)

% maximum number of iterations
maxiter = 1000;

[n,d] = size(X);

X = X - repmat(mean(X),n,1);
% X = X ./ repmat(std(X),n,1);

% calculate jitter weight matrix in D, where multiplication by D is
% equivalent to differencing
D = inv(tril(ones(d)));
A = D;

% A = ones(d,1);
% A(2:2:end) = -1;
% A = A*A';

Ainv = inv(A);


if nargin < 4
    mz = 1;
end


% calculate biased covariance matrix in C, variances in vr
C = X'*X/n;
vr = diag(C);

if nargin < 3
    
    % if no initial conditions are given, choose them at random by scaling
    % standard deviation or variances by a random number in [0,1]

    Z0 = repmat(sqrt(vr),1,mz).*rand(d,mz);
    psi0 = diag(vr.*rand(d,1));
    sigma0 = diag(vr.*rand(d,1));

    %     sigma0 = diag(diag(C));
    %     psi0 = sigma0;
    %     phi0 = diag(rand(1,Qd));


end

% initalized parameter estimates

Z = Z0;
sigma = sigma0;
psi = psi0;

if nargout > 6
    logp = zeros(maxiter,1);
end

iter = 1;

zmat = zeros(maxiter,d,mz);
sigmamat = zeros(maxiter,d);
psimat = sigmamat;

zmat(iter,:,:) = Z;
psimat(iter,:) = diag(psi);
sigmamat(iter,:) = diag(sigma);


% throughout Chat is the estimated covariance matrix, log-likelihood
% follows Bishop (2006)
Chat = psi + A*Z*Z'*A' + D*sigma*D';
logp(iter) = -n/2*(d*log(2*pi) + log(det(Chat)) + trace(inv(Chat)*C));

detpsi = det(psi);
detsigma = det(sigma);

% deltaflg = 1 indicates lack of convergence
deltaflg = 1;

% degflg indicates degeneracy of at least 1 non-global matrix (expected global
% variance = 1 always)
degflg = detpsi > 0 & detsigma > 0;

while iter < maxiter && degflg && deltaflg

    % EM algorithm for estimating parameters as described in Glaze (2008),
    % based in EM for traditional factor analysis in Bishop (2006)
    
    invpsi = inv(psi);

    AZ = A*Z;
    DAZ = [D AZ];
    
    C_uj = [sigma zeros(d,mz);zeros(mz,d) eye(mz)];
    G_uj = inv(inv(C_uj) + DAZ'*invpsi*DAZ);

    uj = X*invpsi*DAZ*G_uj;
    uj2sm = uj'*uj + n*G_uj;

    DAZ = X'*uj * inv(uj2sm);
    Z = Ainv*DAZ(:,d+1:end);
    
    sigma = diag(diag(uj2sm(1:d,1:d))) / n;

    psifull = C + DAZ*uj2sm*DAZ'/n - 2*X'*uj*DAZ'/n;
    psi = diag(diag((psifull)));

    iter = iter + 1;

    psimat(iter,:) = diag(psi);
    zmat(iter,:,:) = Z;
    
    detpsi = det(psi);
    detsigma = det(sigma);
    
    if nargout > 9
        Chat = psi + A*Z*Z'*A';
        logp(iter) = -n/2*(d*log(2*pi) + log(det(Chat)) + trace(inv(Chat)*C));
    end

    
    % settle on estimates if all variance parameters have changed by < 0.01
    % squared-units from the pervious iteration and global weight estimates
    % have changed by < 0.1 units
    
    if max(abs(diff(zmat(iter-1:iter,:)))) < 0.1 && max(abs(diff(psimat(iter-1:iter,:)))) < 0.01 && max(abs(diff(sigmamat(iter-1:iter,:)))) < 0.01
        deltaflg = 0;
    end

    degflg = detpsi > 0 & detsigma > 0;

end

psimat = psimat(1:iter,:);
sigmamat = sigmamat(1:iter,:);
zmat = zmat(1:iter,:);

if nargout > 9
    logp = logp(1:iter);
end

u = uj(:,1:d);
j = uj(:,d+1:end);
