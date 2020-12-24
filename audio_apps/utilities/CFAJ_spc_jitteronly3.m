function [Z,psi,sigma,iter,n_fail,logp,Z0] = CFAJ_spc_jitteronly3(X,bootN,mz)

[n,d] = size(X);

Z = zeros(d,mz,bootN);
psi = zeros(d,bootN);
sigma = psi;
iter = zeros(1,bootN);
logp = iter;

Z0 = Z;
psi0 = Z;
sigma0 = Z;

D = inv(tril(ones(d)));
A = D;

C = cov(X)*n/(n-1);
vr = diag(C);

i = 0;
n_fail = 0;

h = waitbar(0/bootN,'Calculating FA parameters');

while i < bootN

    %     indtmp = ceil(n*rand(n,1));
    
    psi0tmp = diag(vr.*rand(d,1));
    Z0tmp = repmat(sqrt(vr),1,mz).*rand(d,mz);
    sigma0tmp = diag(vr.*rand(d,1));
    
    [Ztmp,psitmp,sigmatmp,itertmp,degflg] = CFA_jitteronly3(X,Z0tmp,psi0tmp,sigma0tmp,mz);

    if itertmp < 1000 && degflg
        i = i + 1;
        
        psi(:,i) = diag(psitmp);
        sigma(:,i) = diag(sigmatmp);
        Z(:,:,i) = Ztmp;
        
        psi0(:,i) = diag(psi0tmp);
        sigma0(:,i) = diag(sigma0tmp);
        Z0(:,:,i) = Z0tmp;

        iter(i) = itertmp;
        
        Chat = A*Ztmp*Ztmp'*A' + psitmp + D*sigmatmp*D';
        logp(i) = -n/2*(d*log(2*pi) + log(det(Chat)) + trace(inv(Chat)*C));
        
        h = waitbar(i/bootN,h,'Calculating FA parameters');
        
    else
        n_fail = n_fail + 1
    end

end

close(h)

