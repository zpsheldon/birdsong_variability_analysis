function [s,c,r] = hiddenmarkov_robgauss_fit(X,A,Pi,mu,sigma,mumu)

[n,d] = size(X);
k = size(A,1);

logdetlambda = zeros(1,k);

a = zeros(n,k);
b = a;
mandist = a;

X1 = repmat(X,[1,1,d]);
X2 = shiftdim(repmat(X',[1,1,d]),1);
Xprod = X1 .* X2;
Xprod = reshape(Xprod,n,d^2);

X2 = reshape(X2,n,d^2);

r = zeros(n,k);

onesvec = ones(n,1);

pidimconst = d*log(2*pi);

for kind = 1:k
    
    lambdatmp = inv(sigma(:,:,kind));
    mulambda = repmat(mu(kind,:)',1,d) .* lambdatmp;
    d1 = Xprod * lambdatmp(:);
    d2 = X2 * mulambda(:);
    mandist(:,kind) = d1 - 2*d2 + trace(mumu(:,:,kind)*lambdatmp);
    
    logdetlambda(kind) = log(det(lambdatmp));
    
    r(:,kind) = max(exp((.5*logdetlambda(kind)-.5*pidimconst)*onesvec - .5*mandist(:,kind)),eps);
    
end

[alphahat,betahat,c,s2sm] = baumwelch(A,r',Pi);
s = (alphahat .* betahat ./ repmat(c,k,1))';
  