function [A_dist,logp_dist] = hiddenmarkov_gauss_latvars_shuff(X,A,mu,sigma,shuffN,initopt)

if nargin < 6
    initopt = 0;
end

A_dist = zeros(shuffN,size(A,1),size(A,2));
logp_dist = zeros(shuffN,1);

N = size(X,1);

if initopt
   initinds = find(X==0); 
end

for shuffind = 1:shuffN
    
    if initopt
        shuffinds = 1:N;
        for initind = 1:length(initinds)-1
           indstmp = shuffinds(initinds(initind)+1:initinds(initind+1)-1);
           shuffindstmp = indstmp(randperm(length(indstmp)));
           shuffinds(indstmp) = shuffindstmp;
        end
    else    
        shuffinds = randperm(N);
    end
    
    Xshuff = X(shuffinds,:);
    [A_dist(shuffind,:,:),gamma,p,logp_dist(shuffind)] = hiddenmarkov_gauss_latvars(Xshuff,A,mu,sigma);
end


