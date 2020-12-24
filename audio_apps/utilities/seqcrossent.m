function crossent = seqcrossent(vc,maxord,words,prl,silind)

t = tabulate(vc);
if nargin == 3
    t = t(t(:,1)~=silind,:);
end

M = t(:,2);
xnm = size(t,1);

P0 = M / sum(M);
minP = 0.0001;

crossent = zeros(1,maxord+1);

pmod = zeros(length(words),1);
for i = 1:length(pmod)
    seq = str2num(words{i});
    pmod(i) = prod(P0(seq))*P0(silind)^2;
end

crossent(1) = -sum(prl .* log2(pmod));

for order = 1:maxord
    
    if nargin == 3
        [transmat,P] = vc2transmat(vc,order,silind);
    else
        [transmat,P] = vc2transmat(vc,order);
    end
    
    X = repmat(sum(transmat,order+1),[ones(1,order) xnm]);
    Pf = transmat ./ X;
    
    dimvc = xnm .^ [0:order-1];
    
    pmod = zeros(length(words),1);
    for i = 1:length(pmod)
        seq = str2num(words{i});
        seq = [silind seq silind];
        
        indmat = cumsum([1:length(seq)-order;ones(order-1,length(seq)-order)],1);
        ind1 = dimvc * seq(indmat);
        ind2 = seq(order+1:end);
        
        pmod(i) = prod(Pf(sub2ind(size(P),ind1,ind2)));
    end
    
    crossent(order+1) = -sum(prl .* log2(pmod));
    
end
