function [Nvc,Pvc,Nvcb,Pvcb] = calcseqstats(seqvc,vc)


seqlen = length(seqvc);

% pad vc with nonsense to avoid searching for values beyond real length
vc = [ones(seqlen,1)*-1;vc(:);ones(seqlen,1)*-1];

vclen = length(vc);

Nvc = zeros(1,seqlen);
Nvcb = Nvc;

inds = find(vc == seqvc(1));
Nvc(1) = length(inds);
for ind = 2:seqlen
    inds = inds(vc(inds+ind-1) == seqvc(ind));
    Nvc(ind) = length(inds);
end

Pvc = Nvc(2:end) ./ Nvc(1:end-1);

if nargout > 2
    
    inds = find(vc(2:end) == seqvc(seqlen))+1;
    Nvcb(seqlen) = length(inds);
    for ind = 2:seqlen
        inds = inds(vc(inds-ind+1) == seqvc(seqlen-ind+1));
        Nvcb(seqlen-ind+1) = length(inds);
    end
    
    Pvcb = Nvcb(1:end-1) ./ Nvcb(2:end);
    
end