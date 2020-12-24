function [scr,w,D,d] = scrdtw(spec1,spec2)

spec1tmp = [zeros(size(spec1,1),3) spec1];
spec2tmp = [zeros(size(spec1,1),3) spec2];

[w,D,d]=dtwDist(spec1tmp,spec2tmp,[1.6 1.6 1]);
scr = prod([size(spec1,1),size(w,1)]) / D(w(1,1),w(1,2));
