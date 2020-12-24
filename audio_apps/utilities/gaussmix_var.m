function [mu,sigma,p,logp,gamma,iter,degflg,logpvc] = gaussmix_var(X,k,mu0,sigma0,p0,crit)

[n,d] = size(X);
C = cov(X);

mu0 = repmat(mean(X),k,1) + (rand(k,d)-.5);
sigma0 = repmat(C,[1,1,k]) + (rand(d,d,k)+1);
p0 = rand(k,1);
p0 = p0 / sum(p0);

B0 = 0.01;
alpha0 = 2;
v0 = n - d - k;
W0 = rand(d);
W0 = W0 * W0'/10 + C;

mu = mu0;
sigma = sigma0;
p = p0;

Winv = zeros(size(sigma));
W = Winv;

maxiter = 10000;
iter = 1;
conv = 0;

if nargout > 7
    logpvc = zeros(1,maxiter);
end

W0inv = inv(W0);

piconst = (2*pi)^(d/2);
invsigma = zeros(size(sigma));
detsigma = zeros(k,1);

expmulambda = zeros(size(X));

logp = 0;

for i = 1:100 %1:20
    i
    for kind = 1:k
        invsigma(:,:,kind) = inv(sigma(:,:,kind));
        sqrtdetsigma(kind) = sqrt(det(sigma(:,:,kind)));
    end
    
    rho = zeros(n,k);
    
    if i == 1
        
        for kind = 1:k
            Xsub = X - repmat(mu(kind,:),n,1);
            
            for nind = 1:n
                rho(nind,kind) = p(kind)*exp(-.5*Xsub(nind,:)*invsigma(:,:,kind)*Xsub(nind,:)')/ (sqrtdetsigma(kind) * piconst);
            end
            
            %         rhotmp = Xsub*invsigma(:,:,kind);
            %         rhotmp = exp(-.5*sum(rhotmp .* Xsub,2));
            %         rho(:,kind) = max(p(kind) * rhotmp / (sqrtdetsigma(kind) * piconst),eps);
            %
        end
        
        
        rhosm = sum(rho,2);
        r = rho ./ repmat(rhosm,1,k);
        
        
    else
        
        for kind = 1:k
            Xsub = X - repmat(m(kind,:),n,1);
            for nind = 1:n
%                 rho(nind,kind)  = exp(explogpi(kind)) * sqrt(exp(exploglambda(kind))) * exp(-d/(2*B(kind)) - .5*v(kind)*Xsub(kind,:)*squeeze(W(:,:,kind))*Xsub(kind,:)');
                rho(nind,kind)  = exp(explogpi(kind)) * sqrt(exp(exploglambda(kind))) * exp(-d/(2*B(kind)) - .5*v(kind)*Xsub(nind,:)*squeeze(W(:,:,kind))*Xsub(nind,:)');
            end
        end
        
        rhosm = sum(rho,2);
        r = rho ./ repmat(rhosm,1,k);
        
    end
    
    N = sum(r);
    xbar = r'*X ./ (repmat(N',1,d));
    
    for kind = 1:k
        sigmatmp = (repmat(sqrt(r(:,kind)),1,d).*(X - repmat(xbar(kind,:),n,1)));
        sigmatmp = sigmatmp'*sigmatmp / N(kind);
        
        if det(sigmatmp) < 0 || isnan(prod(diag(sigmatmp)))
            degflg = 1;
        end
        
        sigma(:,:,kind) = sigmatmp;
    end
    
    B = B0 + N;
    m = (repmat(N',1,d) .* xbar) ./ repmat(B',1,d);
    
    v = v0 + N;
    alpha = alpha0 + N;
    
    exploglambda = zeros(1,k);
    explogpi = exploglambda;
    
    for kind = 1:k
        Winv(:,:,kind) = W0inv + N(kind) * squeeze(sigma(:,:,kind))  + (B0*N(kind)/(B0+N(kind))) * xbar(kind,:)'*xbar(kind,:);
        W(:,:,kind) = inv(Winv(:,:,kind));
              
        exploglambda(kind) = sum(psi((v(kind) + 1 - 1:d) / 2)) + d*log(2) + log(det(W(:,:,kind)));
        explogpi(kind) = psi(alpha(kind)) - psi(sum(alpha));
        
    end
    
 
    %     Winv = shiftdim(repmat(inv(W0),[1,1,k]),2) + repmat(N',d,d) .* shiftdim(sigma,2) + (B0 * N ./ (B0 + N)) .* xbar * xbar';
    
   
    
end
