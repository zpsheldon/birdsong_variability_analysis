function Xshuff = seqshuff(X,initinds)

N = size(X,1);
shuffinds = 1:N;
for initind = 1:length(initinds)-1
    indstmp = shuffinds(initinds(initind)+1:initinds(initind+1)-1);
    shuffindstmp = indstmp(randperm(length(indstmp)));
    shuffinds(indstmp) = shuffindstmp;
end

Xshuff = X(shuffinds,:);