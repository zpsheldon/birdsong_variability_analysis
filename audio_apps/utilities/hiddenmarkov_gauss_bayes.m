function [A,mu,sigma,gamma,Pi,logp,DIC,iter] = hiddenmarkov_gauss_bayes(X,k,alpha,A0,mu0,sigma0)

[n,d] = size(X);
C = cov(X);
piconst = (2*pi)^(d/2);

if nargin < 3
    alpha = 0.001;
end

if nargin < 6 
    
    m = k*rand;
    
    A0 = rand(k);
    A0 = A0 .* repmat(rand(1,k).^(k-m),k,1);
    A0 = A0 ./ repmat(sum(A0,2),1,k);
%     
%     A0 = ones(k)/k;
%     
    mu0 = repmat(mean(X),k,1) .* (rand(k,d)+.5);
    sigma0 = repmat(C,[1,1,k]) .* (rand(d,d,k)+.5);
end
% 
% lambda0 = sigma0;
% mumu0 = sigma0;
% 
% for kind = 1:k
%     lambda0(:,:,kind) = inv(sigma0(:,:,kind));
%     mumu0(:,:,kind) = mu0(kind,:)'*mu0(kind,:);
% end
% 
% mu = mu0;
% lambda = lambda0;
% mumu = mumu0;
% 
% eta0 = d;
% q0 = 0.001;


Pi0 = zeros(1,k);

for kind = 1:k
    Xsub = X(1,:) - mu0(kind,:);
    Pi0(kind) = exp(-.5*diag(Xsub* inv(sigma0(:,:,kind))*Xsub')) / (sqrt(det((sigma0(:,:,kind)))) * piconst);
end

Pi0 = Pi0 / sum(Pi0);

A = A0;
mu = mu0;
sigma = sigma0;
Pi = Pi0;

maxiter = 1000;
iter = 1;
conv = 0;
logp = 0;

invsigma = zeros(d,d,k);
sqrtdetsigma = zeros(1,k);

Betamat = zeros(maxiter,k);

while iter < maxiter && ~conv && ~isnan(A(1,1))
       
    for kind = 1:k
        invsigma(:,:,kind) = inv(max(sigma(:,:,kind),eps));
        sqrtdetsigma(kind) = sqrt(det(max(sigma(:,:,kind),eps)));
    end
    
    B = zeros(n,k);
    
    for kind = 1:k
        Xsub = X - repmat(mu(kind,:),n,1);  
        Btmp = Xsub*invsigma(:,:,kind);
        Btmp = exp(-.5*sum(Btmp .* Xsub,2));        
        B(:,kind) = max(Btmp / (sqrtdetsigma(kind) * piconst),eps);   
    end
    
    [alphahat,betahat,c,gamma2sm] = baumwelch(A,B',Pi);

    gamma = (alphahat .* betahat ./ repmat(c,k,1))';
    
    Betamat(iter,:) = sum(gamma) ./ sum(gamma(:));
    
    Aprev = A;
    
    Atmp = gamma2sm + alpha;
    Asmtmp = repmat(sum(gamma(1:n-1,:))',1,k) + k*alpha;
    
    logA = psi(Atmp) - psi(Asmtmp);
    A = exp(logA);
       
    Pitmp = gamma(1,:)' + alpha;
    Pismtmp = sum(gamma(1,:)') + k*alpha;
    
    logPi = psi(Pitmp) - psi(Pismtmp);
    Pi = exp(logPi);

    k_n = sum(gamma) + k*alpha;
    
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
    
    conv = max(abs(A(:)-Aprev(:))) < 0.0001;
%     conv = (logprev - logp)/logprev < 10e-8;
%     
% 
%     deltap = (logprev - logp)/logprev;
%     conv = deltap < 10e-8;
    
    iter = iter + 1;
    
end

logp = -sum(log(c));

B = zeros(n,k);
for kind = 1:k
    invsigma(:,:,kind) = inv(max(sigma(:,:,kind),eps));
    sqrtdetsigma(kind) = sqrt(det(max(sigma(:,:,kind),eps)));
    
    Xsub = X - repmat(mu(kind,:),n,1);
    Btmp = Xsub*invsigma(:,:,kind);
    Btmp = exp(-.5*sum(Btmp .* Xsub,2));
    B(:,kind) = max(Btmp / (sqrtdetsigma(kind) * piconst),eps);
    
end

Q_A = sum(sum(Atmp .* logA));
Q_Pi = sum(Pitmp .* logPi);
Q_gamma = sum(sum(gamma .* log(B)));

Q = Q_A + Q_Pi + Q_gamma;

DIC = -4*logp + 2*Q;
DIC = [Q_A Q_Pi Q_gamma];

Betamat = Betamat(1:iter-1,:);
