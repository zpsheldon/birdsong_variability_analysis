function inds = findSeq2(x,seq)

%SEARCHES FOR A STANDARD SEQUENCE IN A VECTOR OF INDICES

xlen = length(x);
seqlen = length(seq);

if sign(seq(1))==1
    inds = find(x == seq(1));
else
    inds = find(x ~= -seq(1));
end

inds = inds(inds < xlen);

ind = 2;

while ind <= seqlen && ~isempty(inds)
        
    if sign(seq(ind))==1
        inds = inds(x(inds+ind-1) == seq(ind));
    else
        inds = inds(x(inds+ind-1) ~= -seq(ind));
    end
    
    
    inds = inds(inds < xlen-ind+1);
    ind = ind + 1;
end