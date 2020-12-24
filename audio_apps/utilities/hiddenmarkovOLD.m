function [A,E,p,logp,iter,gamma] = hiddenmarkovOLD(X,k,A0,E0,p0,convcrit)

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
    convcrit = 0.01;
end

A = A0;
E = E0;
p = p0;

maxiter = 1000;
iter = 1;
conv = 0;

logp = 0;


alphahat = zeros(k,n);
betahat = alphahat;

c = zeros(1,n);

Amat = zeros(maxiter,k,k);

while iter < maxiter && ~conv
    
%     alphabar = p .* (X(:,1)'*E)';
%     c(1) = 1 / sum(alphabar);
%     alphahat(:,1) = alphabar * c(1);
%     
%     for nind = 2:n
% %         for kind = 1:k
% %             alphabar(kind) = sum(alphahat(:,nind-1) .* A(:,kind)) * E(seq(nind),kind);
% %         end
%         
%         alphabar = A' * alphahat(:,nind-1) .* (X(:,nind)'*E)';
%         
%         c(nind) = 1 / sum(alphabar);     
%         alphahat(:,nind) = alphabar * c(nind);
%     end
%     
%     betabar = ones(k,1);
%     betahat(:,n) = betabar * c(n);
%     for nind = n-1:-1:1
%         for kind = 1:k
%             betabar(kind) = sum(A(kind,:) .* (X(:,nind+1)'*E) .* betahat(:,nind+1)');
%         end
% 
%         betahat(:,nind) = betabar * c(nind);
%     end
% 
%     
%     gamma2sm = zeros(k);
%    
%     for nind = 1:n-1
%           
%         for kind = 1:k
% %             for kind2 = 1:k
% %                 gamma2sm(nind,kind,kind2) = alphahat(kind,nind) * A(kind,kind2) * (X(:,nind+1)'*E(:,kind2))' * betahat(kind2,nind+1);
% %             end
% %             
%             gamma2sm(kind,:) = gamma2sm(kind,:) + (alphahat(kind,nind) * A(kind,:)' .* (X(:,nind+1)'*E)' .* betahat(:,nind+1))';
%            
%         end
%         
%     end
    
    B = E'*X;
    [alphahat,betahat,c,gamma2sm] = baumwelch(A,B,p);
%     [alphahat,betahat,c,gamma2sm] = baumwelch2(X,A,E,p);
  
    gamma = alphahat .* betahat ./ repmat(c,k,1);
    gamma = gamma';
    
    Aprev = A;
    A = gamma2sm ./ repmat(sum(gamma(1:n-1,:))',1,k);
   
    p = gamma(1,:)';
    p = p / sum(p);
    
    Eprev = E;
    
  	E = X(:,1:n-1)*gamma(1:n-1,:) ./ repmat(sum(gamma(1:n-1,:)),d,1);
    
    
%     logprev = logp;
%     logp = -sum(log(c));
    
    conv = max(abs(A(:)-Aprev(:))) < convcrit && max(abs(E(:)-Eprev(:))) < convcrit;
%     conv = (logprev - logp)/logprev < 10e-8;
    
    iter = iter + 1;
    
end

logp = -sum(log(c));