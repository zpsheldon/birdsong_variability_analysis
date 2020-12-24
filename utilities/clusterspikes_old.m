function [classinds,khat,maxp,varprop,pmean,paramstrct,coeff,scrtot,lat] = clusterspikes_old(X,clustN,eignm,kmax,priork)

[n,d] = size(X);

Xmn = mean(X);
X = X - repmat(Xmn,n,1);
[coeff, scr, lat] = princomp(X);

inds = randperm(n);
inds = inds(1:clustN);

scrtot = scr(inds,1:eignm);

paramstrct = gaussmix_bic(scrtot(:,1:eignm),1:kmax);

paramstrct.bbic = (paramstrct.bic-max(paramstrct.bic)) .* priork;

[dmy,kind] = min(paramstrct.bbic);

khat = paramstrct.k(kind);

vartot = det(cov(scrtot(:,1:eignm)));

pmean = zeros(1,khat);
vark = pmean;

[maxp,classinds] = max(paramstrct.gamma{kind}');

for i = 1:khat
    vark(i) = det(squeeze(paramstrct.sigma{kind}(:,:,i)));
    pmean(i) = mean(maxp(find(classinds==i)));
end

varprop = vark / vartot;


X = scr(:,1:eignm);
mu = paramstrct.mu{kind};
sigma = paramstrct.sigma{kind};
p = paramstrct.p{kind};
n = size(X,1);

gammall = zeros(n,khat);

for kind = 1:khat
    Xsub = X - repmat(mu(kind,:),n,1);
    gammatmp = Xsub*inv(sigma(:,:,kind));
    gammatmp = exp(-.5*sum(gammatmp .* Xsub,2));
    gammall(:,kind) = p(kind) * gammatmp / (sqrt(det(sigma(:,:,kind))) * (2*pi)^(d/2));
    
end

gammasm = sum(gammall,2);
gammall = gammall ./ repmat(gammasm,1,khat);

[maxp,classinds] = max(gammall');