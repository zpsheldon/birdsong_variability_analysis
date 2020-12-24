function [W,psi,sigma,iter,n_fail,logp,W0,psi0,sigma0] = CFAJ_spc_nophi(X,bootN,m_glob)

if nargin < 3
    m_glob = 1;
end

[n,d] = size(X);

W = zeros(d,m_glob,bootN);
sigma = W;
psi = W;
iter = zeros(1,bootN);
logp = iter;

W0 = W;
psi0 = W;
sigma0 = W;

D = inv(tril(ones(d)));
C = cov(X)*n/(n-1);
vr = diag(C);

i = 0;
n_fail = 0;

h = waitbar(0/bootN,'Calculating FA parameters');

while i < bootN

    %     indtmp = ceil(n*rand(n,1));
    
    
    W0tmp = repmat(sqrt(vr),1,m_glob).*rand(d,m_glob);

    sigma0tmp = diag(vr.*rand(d,1));
    psi0tmp = diag(vr.*rand(d,1));

    [Wtmp,sigmatmp,psitmp,itertmp,degflg] = CFA_jitter_nophi(X,m_glob,W0tmp,sigma0tmp,psi0tmp);

    if itertmp < 1000 && degflg
        i = i + 1;

        W(:,:,i) = Wtmp;
        sigma(:,i) = diag(sigmatmp);
        psi(:,i) = diag(psitmp);
        
        W0(:,:,i) = W0tmp;
        sigma0(:,i) = diag(sigma0tmp);
        psi0(:,i) = diag(psi0tmp);
        iter(i) = itertmp;
        
        Chat = Wtmp*Wtmp' + D*sigmatmp*D' + psitmp;
        logp(i) = -n/2*(d*log(2*pi) + log(det(Chat)) + trace(inv(Chat)*C));
        
        h = waitbar(i/bootN,h,'Calculating FA parameters');
        
    else
        n_fail = n_fail + 1
    end

end

close(h)

