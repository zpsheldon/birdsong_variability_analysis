function [A_dist,mu_dist,sigma_dist,logp_dist] = hiddenmarkov_gauss_shuff(X,A,mu,sigma,shuffN,initopt)

if nargin < 6
    initopt = 0;
end

k = size(A,1);
[n,d] = size(X);

A_dist = zeros(shuffN,k,k);
logp_dist = zeros(shuffN,1);
mu_dist = zeros(shuffN,k,d);
sigma_dist = zeros(shuffN,d,d,k);

if initopt
   initinds = find(X==0); 
end

for shuffind = 1:shuffN
    
    if initopt
        shuffinds = 1:n;
        for initind = 1:length(initinds)-1
           indstmp = shuffinds(initinds(initind)+1:initinds(initind+1)-1);
           shuffindstmp = indstmp(randperm(length(indstmp)));
           shuffinds(indstmp) = shuffindstmp;
        end
    else    
        shuffinds = randperm(n);
    end
    
    Xshuff = X(shuffinds,:);
    [Atmp,mutmp,sigmatmp,dmy,dmy,logp_dist(shuffind)] = ...
        hiddenmarkov_gauss(Xshuff,k,A,mu,sigma);
    
    
    A_dist(shuffind,:,:) = Atmp;
    mu_dist(shuffind,:,:) = mutmp;
    sigma_dist(shuffind,:,:,:) = sigmatmp;

end


