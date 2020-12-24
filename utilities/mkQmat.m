function [Q,sequ] = mkQmat(seqarr)

d = length(seqarr);
sequ = unique(seqarr);
du = length(sequ);

Q = zeros(d,du);

indmat = repmat(sequ,d,1);
indmat2 = repmat(seqarr',1,du);

inds = find(strcmp(indmat,indmat2));

Q(inds) = 1;
[i,j] = find(Q);
jtab = tabulate(j);

valinds = jtab(find(jtab(:,2)>1),1);
Q = Q(:,valinds);