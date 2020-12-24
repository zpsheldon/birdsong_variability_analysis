function [W,psi,z,iter,degflg,logp,wmat,psimat] = CFA_Wonly(X,m_glob,W0,psi0)

maxiter = 1000;

if nargin < 2
    m_glob = 1;
end

[n,d] = size(X);

X = X - repmat(mean(X),n,1);
% X = X ./ repmat(std(X),n,1);

C = X'*X/n;
vr = diag(C);

if nargin < 3
    W0 = repmat(sqrt(vr),1,m_glob).*rand(d,m_glob);
    psi0 = diag(vr.*rand(d,1));
end

W = W0;
psi = psi0;

if nargout > 6
    logp = zeros(maxiter,1);
end

iter = 1;

psimat = zeros(maxiter,d);
wmat = zeros(maxiter,d,m_glob);

psimat(iter,:) = diag(psi);
wmat(iter,:,:) = W;


Chat = W*W' + psi;
logp(iter) = -n/2*(d*log(2*pi) + log(det(Chat)) + trace(inv(Chat)*C));

detpsi = det(psi);

deltaflg = 1;
degflg = detpsi > 0;


while iter < maxiter && deltaflg

    iter

    invpsi = diag(1./diag(psi));

    C_z = eye(m_glob);
    G_z = inv(C_z + W'*invpsi*W);

    z = X*invpsi*W*G_z;


    z2sm = z'*z + n*G_z;
    W = (X'*z * inv(z2sm));
    % below involves a sum of inner productss

    psifull = C + W*z2sm*W'/n - 2*X'*z*W'/n;
    psi = diag(diag((psifull)));
    detpsi = det(psi);

    iter = iter + 1;

    wmat(iter,:,:) = W;
    
    if nargout > 4
        Chat = W*W' + psi;
        logp(iter) = -n/2*(d*log(2*pi) + log(det(Chat)) + trace(inv(Chat)*C));
    end

    if max(abs(diff(psimat(iter-1:iter,:)))) < 0.01 && max(max(abs(diff(wmat(iter-1:iter,:,:))))) < 0.01
        deltaflg = 0;
    end

    degflg = detpsi > 0;

end

psimat = psimat(1:iter,:);
wmat = wmat(1:iter,:,:);

if nargout > 6
    logp = logp(1:iter);
end