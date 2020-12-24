function [scr,scr1,scr2,w,D,d] = specdist(spec1,spec2)


pwt = [1.5 1.5 1];
[w,D,d]=dtwDist(spec1,spec2,pwt);

dst = d(sub2ind(size(d),round(w(:,1)),round(w(:,2))));

spec2hat = spec2(:,max(round(sort(w(:,2))),1));

scr = sum(dst) / (size(spec1,1) * size(w,1));
scr1 = corr(spec2hat(:),spec1(:));
scr2 = sqrt(mean((spec1(:)-spec2hat(:)).^2));