function [mu,sigma,p,logp,gamma,iter,degflg,logpvc] = gaussmix(X,k,mu0,sigma0,p0,crit)

[n,d] = size(X);
C = cov(X);

if nargin < 5
    mu0 = repmat(mean(X),k,1) .* (rand(k,d)+1);
    sigma0 = repmat(C,[1,1,k]) .* (rand(d,d,k)+1);
    p0 = rand(k,1);
    p0 = p0 / sum(p0);
end

Xmn = mean(X);
X = X - repmat(Xmn,n,1);

mu = mu0;
sigma = sigma0;
p = p0;

maxiter = 10000;
iter = 1;
conv = 0;

if nargout > 7
    logpvc = zeros(1,maxiter);
end


piconst = (2*pi)^(d/2);
invsigma = zeros(size(sigma));
detsigma = zeros(k,1);

logp = 0;

if nargin < 6
    crit = 0.01;
end

degflg = 0;

while iter < maxiter && ~conv && ~degflg
    
    
    muprev = mu;
    sigmaprev = sigma;
    pprev = p;
    
    for kind = 1:k
        invsigma(:,:,kind) = inv(sigma(:,:,kind));
        sqrtdetsigma(kind) = sqrt(det(sigma(:,:,kind)));
    end
    
    gammasub = zeros(n,k);
    
    for kind = 1:k
        Xsub = X - repmat(mu(kind,:),n,1);
%         for nind = 1:n
%             gammasub(nind,kind) = p(kind)*exp(-.5*Xsub(nind,:)*invsigma(:,:,kind)*Xsub(nind,:)')/ (sqrtdetsigma(kind) * piconst);
%         end
        
        %         gammasub(:,kind) = p(kind) * exp(-.5*diag(Xsub*invsigma(:,:,kind)*Xsub')) / (sqrt(detsigma(kind)) * piconst);
        
        gammatmp = Xsub*invsigma(:,:,kind);
        gammatmp = exp(-.5*sum(gammatmp .* Xsub,2));        
        gammasub(:,kind) = max(p(kind) * gammatmp / (sqrtdetsigma(kind) * piconst),eps);
        
    end
    
    gammasm = sum(gammasub,2);
    gamma = gammasub ./ repmat(gammasm,1,k);
    
    logprev = logp;
    logp = sum(log(gammasm));
    
    if logp > 0
        
        
    end
    
    if nargout > 7
        logpvc(iter) = logp;
    end
    
    k_n = sum(gamma);
    
    mu = gamma'*X ./ (repmat(k_n',1,d));
    for kind = 1:k
        sigmatmp = (repmat(sqrt(gamma(:,kind)),1,d).*(X - repmat(mu(kind,:),n,1)));
        sigmatmp = sigmatmp'*sigmatmp / k_n(kind);
        
        if det(sigmatmp) < 0 || isnan(prod(diag(sigmatmp)))
           degflg = 1; 
        end
            
        sigma(:,:,kind) = sigmatmp;
    end
   
    
    
    p = k_n / n;
    
    %     conv = abs((logprev - logp)/logprev) < 0.001;
    %
    conv = max(abs(p(:)-pprev(:))./pprev(:)) < crit && max(abs(mu(:)-muprev(:)) ./ abs(muprev(:))) < crit && max(abs(sigma(:)-sigmaprev(:)) ./ abs(sigmaprev(:))) < crit;
 
    iter = iter + 1;
    
end

mu = mu + repmat(Xmn,k,1);

if nargout > 7
    logpvc = logpvc(1:iter);
end


