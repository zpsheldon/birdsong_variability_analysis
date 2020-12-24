function [R,featu] = mkRmat(featarr)

d = length(featarr);
featu = unique(featarr);
du = length(featu);

R = zeros(d,du);

indmat = repmat(featu,d,1);
indmat2 = repmat(featarr',1,du);

R(strcmp(indmat,indmat2)) = 1;

D = diff(eye(d));

R = D*R;

[i,j] = find(R);
jtab = tabulate(j);

valinds = jtab(find(jtab(:,2)>2),1);
R = R(:,valinds);