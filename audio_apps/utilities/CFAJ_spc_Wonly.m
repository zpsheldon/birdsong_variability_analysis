function [W,psi,iter,n_fail,logp,W0,psi0] = CFAJ_spc_Wonly(X,bootN,m_glob)

if nargin == 2
    m_glob = 1;
end

[n,d] = size(X);

W = zeros(d,m_glob,bootN);
psi = zeros(d,bootN);
iter = zeros(1,bootN);
logp = iter;

W0 = W;
psi0 = W;

C = cov(X)*n/(n-1);
vr = diag(C);

i = 0;
n_fail = 0;

h = waitbar(0/bootN,'Calculating FA parameters');

while i < bootN

    %     indtmp = ceil(n*rand(n,1));
    
    
    W0tmp = repmat(sqrt(vr),1,m_glob).*rand(d,m_glob);

    psi0tmp = diag(vr.*rand(d,1));
    
    [Wtmp,psitmp,itertmp,degflg] = CFA_Wonly(X,W0tmp,psi0tmp,m_glob);

    if itertmp < 1000 && degflg
        i = i + 1;

        W(:,:,i) = Wtmp;
        psi(:,i) = diag(psitmp);
        
        W0(:,:,i) = W0tmp;
        psi0(:,i) = diag(psi0tmp);
    
        iter(i) = itertmp;
        
        Chat = Wtmp*Wtmp' + psitmp;
        logp(i) = -n/2*(d*log(2*pi) + log(det(Chat)) + trace(inv(Chat)*C));
        
        h = waitbar(i/bootN,h,'Calculating FA parameters');
        
    else
        n_fail = n_fail + 1
    end

end

close(h)

