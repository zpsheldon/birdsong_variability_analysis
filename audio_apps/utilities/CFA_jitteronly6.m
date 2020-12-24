function [sigma,psi,iter,degflg,u,logp,sigmamat,psimat] = CFA_jitteronly6(X,sigma0,psi0)

maxiter = 1000;

[n,d] = size(X);

X = X - repmat(mean(X),n,1);
% X = X ./ repmat(std(X),n,1);

D = -diff(eye(d))';

C = X'*X/n;
vr = diag(C);

if nargin < 4
    sigma0 = diag([vr(1:d-1)+vr(2:d)].*rand(d-1,1));
    psi0 = diag(vr.*rand(d,1));
end

sigma = sigma0;
psi = psi0;

if nargout > 5
    logp = zeros(maxiter,1);
end

iter = 1;

sigmamat = zeros(maxiter,d-1);
psimat = zeros(maxiter,d);

sigmamat(iter,:) = diag(sigma);
psimat(iter,:) = diag(psi);

Chat = D*sigma*D' + psi;
logp(iter) = -n/2*(d*log(2*pi) + log(det(Chat)) + trace(inv(Chat)*C));

detsigma = det(sigma);
detpsi = det(psi);

deltaflg = 1;
degflg = detsigma > 0 & detpsi > 0;


while iter < maxiter && degflg && deltaflg

    %     iter

    invpsi = inv(psi);
    
    C_zu = sigma;
    G_zu = inv(inv(C_zu) + D'*invpsi*D);

    u = X*invpsi*D*G_zu;

    u2sm = u'*u + n*G_zu;

    sigma = diag(diag(u2sm)) / n;


    % below involves a sum of inner productss

    psifull = C + D*u2sm*D'/n - 2*X'*u*D'/n;
    psi = diag(diag((psifull)));

    detsigma = det(sigma);
    detpsi = det(psi);

    iter = iter + 1;

    sigmamat(iter,:) = diag(sigma);
    psimat(iter,:) = diag(psi);
    
    if nargout > 5
        Chat = D*sigma*D' + psi;
        logp(iter) = -n/2*(d*log(2*pi) + log(det(Chat)) + trace(inv(Chat)*C));
    end

    if max(abs(diff(sigmamat(iter-1:iter,:)))) < 0.001 & max(abs(diff(psimat(iter-1:iter,:)))) < 0.001
        deltaflg = 0;
    end

    degflg = detsigma > 0 & detpsi > 0;

end

psimat = psimat(1:iter,:);
sigmamat = sigmamat(1:iter,:);

if nargout > 6
    logp = logp(1:iter);
end


