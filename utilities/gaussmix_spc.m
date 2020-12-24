function [mu,sigma,p,gamma,logp,iter,n_fail] = gaussmix_spc(X,k,initN)

if nargin < 3
    initN = 100;
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
logprev = -inf;

itot = 0;

while i < initN & itot < 5*initN
    
    %     indtmp = ceil(n*rand(n,1));
    
    
    mu0tmp = stdXmat .* (rand(k,d)-.5);
    
%     sigmacholtmp = Smat .* (.5 + rand(d,d,k));
%     
%     sigma0tmp = zeros(d,d,k);
%     
%     for kind = 1:k
%         sigmachol = triu(squeeze(sigmacholtmp(:,:,kind)));
%         sigma0tmp(:,:,kind) = sigmachol * sigmachol';
%         
%         if det(squeeze(sigma0tmp(:,:,kind))) <= 0
%             
%         end
%         
%     end
%     
    
    sigma0tmp = Smat;
    
    p0tmp = rand(1,k);
    p0tmp = p0tmp / sum(p0tmp);
    
    [mutmp,sigmatmp,ptmp,logptmp,gammatmp,itertmp,degflg] = gaussmix(X,k,mu0tmp,sigma0tmp,p0tmp,0.001);
    
    itot = itot+1;
    
    if itertmp < 10000 && ~degflg
        
        i = i + 1;
        h = waitbar(i/initN,h,'Calculating parameters');
        
        if logptmp > logprev
            mu = mutmp;
            sigma = sigmatmp;
            p = ptmp;
            logp = logptmp;
            iter = itertmp;
            gamma = gammatmp;
        end
        
        logprev = logptmp;
        
    else
        n_fail = n_fail + 1
    end
    
end

if i == 0
    logp = -inf;
    mu = [];
    sigma = [];
    p = [];
    iter = 0;
    gamma = [];
end

close(h)