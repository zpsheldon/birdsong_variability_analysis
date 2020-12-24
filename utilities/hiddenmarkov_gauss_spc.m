function [A,Pi,mu,sigma,gamma,logp,BIC,m,logpvc,logpall] = hiddenmarkov_gauss_spc(X,k,initN,zeroinds)

if nargin == 2
    initN = 10;
end

if nargin < 4
    zeroinds = [];
end

i = 0;
logp = -inf;
BIC = inf;

[n,d] = size(X);

logpall = zeros(initN,1);

while i < initN

    [Atmp,mutmp,sigmatmp,gammatmp,Pitmp,logptmp,iter,logpvctmp] = hiddenmarkov_gauss(X,k,zeroinds);
    ptmp = max(sum(gammatmp) / sum(gammatmp(:)),eps);
    mhat = prod(ptmp .^ -ptmp);
    
    khat = mhat*(mhat-1) + mhat*d + mhat*d*(d+1)/2 + mhat - 1;
    
    BICtmp = -2*logptmp + log(n)*khat;


    if i==1 || logptmp > logp
        
%     if i==1 || BICtmp < BIC
       A = Atmp;
       mu = mutmp;
       sigma = sigmatmp;
       gamma = gammatmp;
       logp = logptmp;
       Pi = Pitmp;
       BIC = BICtmp;
       m = mhat;
       logpvc = logpvctmp;
    end
    
    i = i + 1;
        
    logpall(i) = logptmp;
end

 [A,mu,sigma,gamma,Pi,logp,iter,logpvc] = hiddenmarkov_gauss(X,k,zeroinds,A,Pi,mu,sigma,1);
 