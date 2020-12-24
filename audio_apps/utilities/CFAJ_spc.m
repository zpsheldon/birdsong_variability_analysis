function [W,psi,phi,sigma,iter,n_fail,logp,W0,psi0,sigma0,phi0] = CFAJ_spc(X,Q,bootN,m_glob)

if nargin < 4
    m_glob = 1;
end

[n,d] = size(X);
Qd = size(Q,2);

W = zeros(d,m_glob,bootN);
sigma = zeros(d,bootN);
psi = sigma;
phi = zeros(Qd,bootN);
iter = zeros(1,bootN);
logp = iter;

W0 = W;
psi0 = sigma;
sigma0 = sigma;
phi0 = phi;

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

    [rows,cols] = find(Q);
    rowsu = unique(rows);

    phi0tmp = diag(vr(rowsu(1:Qd)) .* rand(Qd,1));
    
    [Wtmp,sigmatmp,psitmp,phitmp,itertmp,degflg,z,u,q] = CFA_jitter(X,Q,m_glob,W0tmp,sigma0tmp,psi0tmp,phi0tmp);

    if itertmp < 1000 && degflg
        i = i + 1;

        W(:,:,i) = Wtmp;
        sigma(:,i) = diag(sigmatmp);
        psi(:,i) = diag(psitmp);
        
        W0(:,:,i) = W0tmp;
        sigma0(:,i) = diag(sigma0tmp);
        psi0(:,i) = diag(psi0tmp);
    
        if Qd > 0
            phi(:,i) = diag(phitmp);
            phi0(:,i) = diag(phi0tmp);
        end

        iter(i) = itertmp;
        
        Chat = Wtmp*Wtmp' + D*sigmatmp*D' + Q*phitmp*Q' + psitmp;
        logp(i) = -n/2*(d*log(2*pi) + log(det(Chat)) + trace(inv(Chat)*C));
        
        h = waitbar(i/bootN,h,'Calculating FA parameters');
        
    else
        n_fail = n_fail + 1
    end

end

close(h)

