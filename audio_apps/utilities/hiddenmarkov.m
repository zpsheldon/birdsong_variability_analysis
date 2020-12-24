function [A,E,p,logp,iter,gamma] = hiddenmarkov(X,k,A0,E0,p0,convcrit)

[d,n] = size(X);

if nargin < 5
    
    A0 = rand(k);
    A0 = A0 ./ repmat(sum(A0,2),1,k);
    
    E0 = rand(d,k);
    E0 = E0 ./ repmat(sum(E0,1),d,1);
    
    p0 = rand(k,1);
    p0 = p0 / sum(p0);
    
end

if nargin < 6
    convcrit = 0.001;
end

A = A0;
E = E0;
p = p0;

maxiter = 1000;
iter = 1;
conv = 0;
logp = 0;

while iter < maxiter && ~conv
       
    B = E'*X;
    [alphahat,betahat,c,gamma2sm] = baumwelch(A,B,p);

    gamma = (alphahat .* betahat ./ repmat(c,k,1))';
    
    Aprev = A;
    A = gamma2sm ./ repmat(sum(gamma(1:n-1,:))',1,k);
   
    p = gamma(1,:)';
    p = p / sum(p);
    
    Eprev = E;
    
  	E = X(:,1:n-1)*gamma(1:n-1,:) ./ repmat(sum(gamma(1:n-1,:)),d,1);
    
    
    logprev = logp;
    logp = -sum(log(c));
    
% %     conv = max(abs(A(:)-Aprev(:))) < convcrit && max(abs(E(:)-Eprev(:))) < convcrit;
%     conv = (logprev - logp)/logprev < 10e-8;
%     

    deltap = (logprev - logp)/logprev
    conv = deltap < 0.001;
    
    iter = iter + 1
    
end

logp = -sum(log(c));