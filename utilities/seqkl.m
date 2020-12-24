function klvc = seqkl(vc,maxord,silind)

t = tabulate(vc);
if nargin == 3
    t = t(t(:,1)~=silind,:);
end

M = t(:,2);
xnm = size(t,1);

Pk = M / sum(M);
minP = 0.0001;

klvc = zeros(1,maxord);

for order = 1:maxord
    
    if nargin == 3
        [transmat,P] = vc2transmat(vc,order,silind);
    else
        [transmat,P] = vc2transmat(vc,order);
    end
    
    if order == 1
        P0 = Pk * Pk';
        P1 = P;      
    else
        P1tmp = repmat(Pprev,[ones(1,order) xnm]);
        P2tmp = shiftdim(repmat(Pprev2,[ones(1,order) xnm]),order);
        P0 = P1tmp .* P2tmp;
    end
    
    X = repmat(sum(transmat,order+1),[ones(1,order) xnm]);
    Pprev2 = transmat ./ X;
    Pprev = P;
    
    kl = max(P,minP) .* log2(max(P,minP) ./ max(P0,minP));
    klvc(order) = sum(kl(:));
end
