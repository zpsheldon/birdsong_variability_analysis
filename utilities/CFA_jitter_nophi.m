function [W,sigma,psi,iter,degflg,z,u,logp,wmat,sigmamat,psimat] = CFA_jitter_nophi(X,m_glob,W0,sigma0,psi0)

maxiter = 1000;

[n,d] = size(X);

X = X - repmat(mean(X),n,1);
% X = X ./ repmat(std(X),n,1);

D = inv(tril(ones(d)));

C = X'*X/n;
vr = diag(C);

if nargin < 4

    W0 = repmat(sqrt(vr),1,m_glob).*rand(d,m_glob);

    sigma0 = diag(vr.*rand(d,1));
    psi0 = diag(vr.*rand(d,1));
    
end

W = W0;
sigma = sigma0;
psi = psi0;

if nargout > 6
    logp = zeros(maxiter,1);
end

iter = 1;

sigmamat = zeros(maxiter,d);
psimat = sigmamat;
wmat = zeros(maxiter,d,m_glob);

sigmamat(iter,:) = diag(sigma);
psimat(iter,:) = diag(psi);
wmat(iter,:,:) = W;


Chat = W*W' + D*sigma*D' + psi;
logp(iter) = -n/2*(d*log(2*pi) + log(det(Chat)) + trace(inv(Chat)*C));

detsigma = det(sigma);
detpsi = det(psi);

deltaflg = 1;
degflg = detsigma > 0 & detpsi > 0;


while iter < maxiter && degflg && deltaflg

    %     iter

    invpsi = inv(psi);
    WD = [W D];

    C_zu = [eye(m_glob) zeros(m_glob,d);zeros(d,m_glob) sigma];
    G_zu = inv(inv(C_zu) + WD'*invpsi*WD);

    zu = X*invpsi*WD*G_zu;


    zu2sm = zu'*zu + n*G_zu;

    %     sigma = diag(diag(u2sm)) / n;

    sigma = diag(diag(zu2sm(m_glob+1:m_glob+d,m_glob+1:m_glob+d))) / n;

    WD = (X'*zu * inv(zu2sm));
    
    W = WD(:,1:m_glob);

    % below involves a sum of inner productss

    psifull = C + WD*zu2sm*WD'/n - 2*X'*zu*WD'/n;
    psi = diag(diag((psifull)));



    detsigma = det(sigma);
    detpsi = det(psi);

    iter = iter + 1;

    sigmamat(iter,:) = diag(sigma);
    psimat(iter,:) = diag(psi);
    wmat(iter,:,:) = W;
    
    if nargout > 7
        Chat = W*W' + D*sigma*D' + psi;
        logp(iter) = -n/2*(d*log(2*pi) + log(det(Chat)) + trace(inv(Chat)*C));
    end

    if max(abs(diff(sigmamat(iter-1:iter,:)))) < 0.01 & max(abs(diff(psimat(iter-1:iter,:)))) < 0.01 ...
            &  max(max(abs(diff(wmat(iter-1:iter,:,:))))) < 0.1 
        deltaflg = 0;
    end

    degflg = detsigma > 0 & detpsi > 0;

end

psimat = psimat(1:iter,:);
sigmamat = sigmamat(1:iter,:);
wmat = wmat(1:iter,:,:);

if nargout > 6
    logp = logp(1:iter);
end

z = zu(:,1:m_glob);
u = zu(:,m_glob+1:m_glob+d);

