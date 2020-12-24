function [mu,sigma,p,v,u,L,W,s,iter,n_fail,mumat,sigmamat,pmat,vmat,Lvc] = robgaussmix_spc(X,k,initN)

if nargin < 3
    initN = 30;
end

[n,d] = size(X);

S = cov(X);
% Smat = repmat(chol(S),[1,1,k]);
Smat = repmat(S,[1,1,k]);
muX = mean(X);
muXmat = repmat(muX,[k,1]);
stdX = std(X);
stdXmat = repmat(stdX,[k,1]);

i = 0;
n_fail = 0;

h = waitbar(0/initN,'Calculating parameters');
Lmax = -inf;

itot = 0;

v0tmp = 0.1*ones(1,k);

if nargout >= 11
    mumat = zeros(k,d,initN);
    sigmamat = zeros(d,d,k,initN);
    pmat = zeros(k,initN);
    vmat = zeros(k,initN);
    Lvc = zeros(1,initN);
end

while i < initN && itot < 5*initN
 
    mu0tmp = stdXmat .* (rand(k,d)-.5);
    sigma0tmp = Smat;
    
    p0tmp = rand(1,k);
    p0tmp = p0tmp / sum(p0tmp);
    
    [mutmp,sigmatmp,ptmp,vtmp,Ltmp,Wtmp,itertmp,s,u,mumu,logu,logdetlambda] = robgaussmix(X,k,mu0tmp,sigma0tmp,p0tmp,v0tmp,1);
    
    itot = itot+1;

    if itertmp < 1000
        
        i = i + 1;
        
        
            
        if nargout >= 11
            mumat(:,:,i) = mutmp;
            sigmamat(:,:,:,i) = sigmatmp;
            pmat(:,i) = ptmp;
            vmat(:,i) = vtmp;
            Lvc(i) = Ltmp;
        end
        
        
        h = waitbar(i/initN,h,'Calculating parameters');
        
        if Ltmp > Lmax
            mu0 = mu0tmp;
            sigma0 = sigma0tmp;
            p0 = p0tmp;
            Lmax = Ltmp;           
        end
        
    else
        n_fail = n_fail + 1
    end
    
end

if i == 0
    L = -inf;
    mu = [];
    sigma = [];
    p = [];
    iter = 0;
    s = [];
    v = [];
    u = [];
    W = [];
else
   [mu,sigma,p,v,L,W,iter,s,u] = robgaussmix(X,k,mu0,sigma0,p0,v0tmp,0); 
end

close(h)