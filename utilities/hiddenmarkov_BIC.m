function [A,E,p,logpvc,bicvc] = hiddenmarkov_BIC(X,kvc)

if nargin == 1
    kvc = 1:5;
end

knm = length(kvc);
logpvc = zeros(1,knm);
mvc = logpvc;

[d,n] = size(X);

bicflg = 1;

for kind = 1:knm
    kind
    [A{kind},E{kind},p{kind},dmy,dmy,logpvc(kind)] = hiddenmarkov_spc(X,kvc(kind),100);
%     mvc(kind) = (kvc(kind) + d + 1) * (kvc(kind)-1);
    mvc(kind) = kvc(kind)*(kvc(kind)-1) - 1 + (kvc(kind)-1)*(d-1);
    bicvc(kind) = -2 * logpvc(kind) + log(n) * mvc(kind);
end



% [bic,minind] = min(bicvc);
% A = A{minind};
% E = E{minind};
% p = p{minind};
% logp = logpvc(minind);