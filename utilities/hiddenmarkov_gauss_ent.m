function [A,mu,sigma,gamma] = hiddenmarkov_gauss_ent(X,kvc)

entvc = zeros(1,length(kvc));
mvc = entvc;
bicvc = entvc;
logpvc = entvc;

[n,d] = size(X);

for kind = 1:length(kvc)
    [A{kind},mu{kind},sigma{kind},gamma{kind},logpvc(kind)] = hiddenmarkov_gauss_spc(X,kvc(kind),10);

    gammatmp = gamma{kind};
    Atmp = A{kind};
    
    entvc(kind) = -sum(sum(Atmp' .* log(max(Atmp',.00001))) .* sum(gammatmp))/sum(gammatmp(:));
    
    
    mvc(kind) = kvc(kind)*(kvc(kind)-1) - 1 + (kvc(kind)-1)*(d-1);
    bicvc(kind) = -2 * logpvc(kind) + log(n) * mvc(kind);
    
    kvc(kind)
    
end



[bic,minind] = min(entvc);
% A = A{minind};
% E = E{minind};
% p = p{minind};
% logp = logpvc(minind);


A = A{minind};
mu = mu{minind};
sigma = sigma{minind};
gamma = gamma{minind};