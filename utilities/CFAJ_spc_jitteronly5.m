function [psi,sigma,iter,n_fail,logp,psi0,sigma0] = CFAJ_spc_jitteronly5(X,bootN)

[n,d] = size(X);

sigma = zeros(d+1,bootN);
psi = zeros(d,bootN);
iter = zeros(1,bootN);
logp = iter;

psi0 = psi;
sigma0 = sigma;

D = diff(eye(d+1));
C = cov(X)*n/(n-1);
vr = diag(C);

i = 0;
n_fail = 0;

h = waitbar(0/bootN,'Calculating FA parameters');

while i < bootN

    %     indtmp = ceil(n*rand(n,1));
    
    
    sigma0tmp = diag([vr(1);vr(1:d-1)+vr(2:d);vr(d)].*rand(d+1,1));
    psi0tmp = diag(vr.*rand(d,1));

    [sigmatmp,psitmp,itertmp,degflg] = CFA_jitteronly5(X,sigma0tmp,psi0tmp);

    if itertmp < 1000 && degflg
        i = i + 1;

        sigma(:,i) = diag(sigmatmp);
        psi(:,i) = diag(psitmp);
        
        sigma0(:,i) = diag(sigma0tmp);
        psi0(:,i) = diag(psi0tmp);
        iter(i) = itertmp;
        
        Chat = D*sigmatmp*D' + psitmp;
        logp(i) = -n/2*(d*log(2*pi) + log(det(Chat)) + trace(inv(Chat)*C));
        
        h = waitbar(i/bootN,h,'Calculating FA parameters');
        
    else
        n_fail = n_fail + 1
    end

end

close(h)

