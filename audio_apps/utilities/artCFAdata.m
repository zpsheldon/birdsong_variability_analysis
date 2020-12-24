function X = artCFAdata(W,psi,sigma,N)

jitter = sqrt(sigma);
ind = sqrt(psi);

M = length(W(:));

jtmp = zeros(N,M-1);

for h = 1:M-1
    jtmp(:,h) = normrnd(0,jitter(h),N,1);
end

j2 = [normrnd(0,mean(jitter),N,1) jtmp normrnd(0,mean(jitter),N,1)];
j2 = diff(j2')';

n = zeros(N,M);
for h = 1:M 
    n(:,h) = normrnd(0,ind(h),N,1);
end

g = normrnd(0,1,N,1)*W(:)';

X = n + j2 + g;