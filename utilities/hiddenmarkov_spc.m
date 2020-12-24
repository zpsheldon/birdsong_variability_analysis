function [A,E,p,iter,n_fail,logp] = hiddenmarkov_spc(X,k,initN,lropt)

if nargin == 2
    initN = 100;
end

if nargin < 4
    lropt = 0;
end

[d,n] = size(X);

A = zeros(k,k,initN);
E = zeros(d,k,initN);
p = zeros(k,initN);
iter = zeros(1,initN);
logp = iter;

i = 0;
n_fail = 0;

% h = waitbar(0/initN,'Calculating FA parameters');

while i < initN
    
    %     indtmp = ceil(n*rand(n,1));
        
%     A0tmp = 2*eye(k);
%     A0tmp(find(~A0tmp)) = rand(1,k^2 - k);
    
%     A0tmp = rand(k);
%     A0tmp = A0tmp ./ repmat(sum(A0tmp,2),1,k);
%     
%     E0tmp = rand(d,k);
%     E0tmp = E0tmp ./ repmat(sum(E0tmp,1),d,1);
%     
%     p0tmp = rand(k,1);
%     p0tmp = p0tmp / sum(p0tmp);
%     
%     
    A0tmp = rand(k);
    A0tmp(1) = 0;
%     A0tmp(2:end,2:end) = triu(A0tmp(2:end,2:end));
    A0tmp = A0tmp ./ repmat(sum(A0tmp,2),1,k);
    
    E0tmp = rand(d,k);
    E0tmp(:,1) = [zeros(d-1,1);1];
    E0tmp(end,2:end) = 0;
    E0tmp = E0tmp ./ repmat(sum(E0tmp,1),d,1);
    
    p0tmp = [1;zeros(k-1,1)];

    [Atmp,Etmp,ptmp,logptmp,itertmp] = hiddenmarkov(X,k,A0tmp,E0tmp,p0tmp);
   
    if itertmp < 1000
        i = i + 1
        
        A(:,:,i) = Atmp;
        E(:,:,i) = Etmp;
        p(:,i) = ptmp;
        iter(i) = itertmp;    
        logp(i) = logptmp;
        
        %         h = waitbar(i/initN,h,'Calculating FA parameters');
        
    else
        n_fail = n_fail + 1
    end
    
end

[logp,maxind] = max(logp);
A = A(:,:,maxind);
E = E(:,:,maxind);
p = p(:,maxind);
iter = iter(maxind);
% 
% alpha = zeros(k,n);
% beta = alpha;
% 
% alphasm = zeros(1,n);
% betasm = alphasm;
% 
% [seq,dmy] = find(X);
% 
% alpha(:,1) = p' .* E(seq(1),:);
% alphasm(1) = sum(alpha(:,1));
% alpha(:,1) = alpha(:,1) / alphasm(1);
% for nind = 2:n
%     alpha(:,nind) = alpha(:,nind-1)'*A .* E(seq(nind),:);
%     alphasm(nind) = sum(alpha(:,nind));
%     alpha(:,nind) = alpha(:,nind) / alphasm(nind);
% end
% 
% %     alphas = alpha.*repmat(cumprod(alphasm),size(alpha,1),1);
% 
% beta(:,end) = ones(k,1);
% betasm(end) = sum(beta(:,end));
% beta(:,end) = beta(:,end) / betasm(end);
% for nind = n-1:-1:1
%     beta(:,nind) = beta(:,nind+1)'*A .* E(seq(nind+1),:);
%     betasm(nind) = sum(beta(:,nind));
%     beta(:,nind) = beta(:,nind) / betasm(nind);
% end
% 
% %     betas = beta.*repmat(fliplr(cumprod(fliplr(betasm))),size(beta,1),1);
% 
% gamma = alpha .* beta;
% gamma = gamma ./ repmat(sum(gamma),k,1);
% 

% close(h)
