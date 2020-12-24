function [transmat,P,xu,kl] = vc2transmat(vc,order,silind)

vc = vc(:)';

if nargin == 1
    order = 1;
end

t = tabulate(vc);
xu = t(:,1);
M = t(:,2);
xnm = length(xu);

transmat = zeros(xnm*ones(1,order+1));

dimtrans = xnm .^ [0:order];
indtrans = cumsum([1:(length(vc)-order);ones(order,length(vc)-order)]);
vclen = length(vc)-order;

prinds = dimtrans * [vc(indtrans) - [zeros(1,vclen);ones(order,vclen)]];

t = tabulate(prinds);
itmp2 = find(t(:,2));
transmat(t(itmp2,1)) = transmat(t(itmp2,1)) + t(itmp2,2);

if nargin == 3
   indkern = setdiff(1:length(M),silind);
   exptmp = 'indkern,';
   exp = [];
   for i = 1:order+1
      exp = [exp exptmp];
   end
   
   exp = exp(1:end-1);
   
   eval(['transmat = transmat(' exp ');']);
   
end

P = transmat / sum(transmat(:));

if nargout==4
    Ptmp = repmat(M/sum(M),[1 xnm*ones(1,order)]);
    P0 = Ptmp;
    for i = 1:order
        P0 = P0 .* shiftdim(Ptmp,i);
    end
    
    minP = 0.00001;
    
    kl = max(P,minP) .* log2(max(P,minP) ./ P0);
    kl = sum(kl(:));
end
