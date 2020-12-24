function [A,Pi,mu,sigma,v,s,L,BIC,m,u,logu,mumu] = hiddenmarkov_robgauss_spc(X,k,initN,alphak,zeroinds)

if nargin == 2
    initN = 10;
end

if nargin < 4
    alphak = 0.001;
end

if nargin < 5
    zeroinds = [];
end

i = 0;
L = -inf;

BIC = inf;

[n,d] = size(X);

while i < initN
    
    [Atmp,mutmp,sigmatmp,Pitmp,vtmp,Ltmp,DICtmp,Wtmp,dmy,stmp,LLvc,utmp,mumutmp,logutmp] = hiddenmarkov_robgauss(X,k,alphak,zeroinds,1);

    eta0 = sum(stmp) + k*alphak;
    etasm0 = sum(eta0);
    logBeta0 = psi(eta0) - psi(etasm0);
    
    Beta0 = exp(logBeta0);
    
    mhat = exp(-sum(Beta0 .* logBeta0));
    khat = mhat*(mhat-1) + mhat*d + mhat*d*(d+1)/2 + mhat;
    
    BICtmp = -2*Ltmp + log(n)*khat;
    
    Ltmp
    
%     if i==1 || BICtmp < BIC
   if i==1 || Ltmp > L
        A = Atmp;
        mu = mutmp;
        Pi = Pitmp;
        sigma = sigmatmp;
        v = vtmp;
        s = stmp;
        L = Ltmp;
        m = mhat;
        u = utmp;
        mumu = mumutmp;
        logu = logutmp;
        BIC = BICtmp;
    end
    
    i = i + 1
end

[A,mu,sigma,Pi,v,L,DIC,W,dmy,s,LLvc2,u,mumu,logu] = hiddenmarkov_robgauss(X,k,alphak,zeroinds,0,A,Pi,mu,sigma,v,u,logu,mumu);