function paramstrct = gaussmix_bic(X,kvc)

if nargin < 2
    kvc = 1:5;
end

paramstrct.k = kvc;

kind = 1;

[n,d] = size(X);

for k = kvc
    
    k
    
    [mu,sigma,p,gamma,logp,iter,n_fail] = gaussmix_spc(X,k,10);
    
    paramstrct.mu{kind} = mu;
    paramstrct.sigma{kind} = sigma;
    paramstrct.p{kind} = p;
    paramstrct.logp(kind) = logp;
    paramstrct.gamma{kind} = gamma;
    
    paramstrct.paramnm(kind) = kvc(kind)*d + kvc(kind)*d*(d+1)/2 + kvc(kind) - 1;
    
    kind = kind + 1;
    
end

paramstrct.bic = -2 * paramstrct.logp + paramstrct.paramnm * log(n);