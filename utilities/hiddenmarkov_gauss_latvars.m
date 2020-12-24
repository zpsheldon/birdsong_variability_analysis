function [A_X,gamma,p_X,logp] = hiddenmarkov_gauss_latvars(X,A,mu,sigma)

[n,d] = size(X);
k = size(A,1);
ptmp = zeros(1,k);
piconst = (2*pi)^(d/2);

for kind = 1:k
    Xsub = X(1,:) - mu(kind,:);
    ptmp(kind) = max(exp(-.5*diag(Xsub* inv(sigma(:,:,kind))*Xsub')) / (sqrt(det((sigma(:,:,kind)))) * piconst),eps);
end
p = ptmp / sum(ptmp);

B = zeros(n,k);

for kind = 1:k
    Xsub = X - repmat(mu(kind,:),n,1);
    Btmp = Xsub*inv(sigma(:,:,kind));
    Btmp = exp(-.5*sum(Btmp .* Xsub,2));
    B(:,kind) = max(Btmp / (sqrt(det((sigma(:,:,kind)))) * piconst),eps);
end

[alphahat,betahat,c,gamma2sm] = baumwelch(A,B',p);

gamma = (alphahat .* betahat ./ repmat(c,k,1))';

A_X = gamma2sm ./ repmat(sum(gamma(1:n-1,:))',1,k);

p_X = gamma(1,:)';
p_X = p_X / sum(p_X);

logp = -sum(log(c));