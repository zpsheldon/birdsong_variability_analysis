function [A,mu,sigma,gamma,p,logp,iter] = hiddenmarkov_gauss_alt(X,k,A0,E0,p0,convcrit)

[n,d] = size(X);
C = cov(X);

if nargin < 5
    
    A0 = rand(k);
    A0 = A0 ./ repmat(sum(A0,2),1,k);
    
    mu0 = repmat(mean(X),k,1) .* (rand(k,d)+1);
    sigma0 = repmat(C,[1,1,k]) .* (rand(d,d,k)+1);
        
    p0 = rand(k,1);
    p0 = p0 / sum(p0);
    
end

if nargin < 6
    convcrit = 0.001;
end

A = A0;
mu = mu0;
sigma = sigma0;
p = p0;

maxiter = 1000;
iter = 1;
conv = 0;
logp = 0;

invsigma = zeros(d,d,k);
sqrtdetsigma = zeros(1,k);

piconst = (2*pi)^(d/2);

while iter < maxiter && ~conv
       
%     B = mu'*X;
    
   
    for kind = 1:k
        invsigma(:,:,kind) = inv(sigma(:,:,kind));
        sqrtdetsigma(kind) = sqrt(det(sigma(:,:,kind)));
    end
    
    B = zeros(n,k);
    
    for kind = 1:k
        Xsub = X - repmat(mu(kind,:),n,1);
%         for nind = 1:n
%             gammasub(nind,kind) = p(kind)*exp(-.5*Xsub(nind,:)*invsigma(:,:,kind)*Xsub(nind,:)')/ (sqrtdetsigma(kind) * piconst);
%         end
        
        %         gammasub(:,kind) = p(kind) * exp(-.5*diag(Xsub*invsigma(:,:,kind)*Xsub')) / (sqrt(detsigma(kind)) * piconst);
        
        Btmp = Xsub*invsigma(:,:,kind);
        Btmp = exp(-.5*sum(Btmp .* Xsub,2));        
        B(:,kind) = max(Btmp / (sqrtdetsigma(kind) * piconst),eps);
        
    end
    
    [alphahat,betahat,c,gamma2sm] = baumwelch(A,B',p);

    gamma = (alphahat .* betahat ./ repmat(c,k,1))';
    
    Aprev = A;
    A = gamma2sm ./ repmat(sum(gamma(1:n-1,:))',1,k);
   
    p = gamma(1,:)';
    p = p / sum(p);

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
    
    
    logprev = logp;
    logp = -sum(log(c));
    
% %     conv = max(abs(A(:)-Aprev(:))) < convcrit && max(abs(E(:)-Eprev(:))) < convcrit;
%     conv = (logprev - logp)/logprev < 10e-8;
%     

    deltap = (logprev - logp)/logprev;
    conv = deltap < 10e-8;
    
    iter = iter + 1
    
end

logp = -sum(log(c));