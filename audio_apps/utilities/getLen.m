function [l,onset,offset] = getLen(s)

winlen = 256;
winadv = 32;

% winadv = 128;

s = s(:)';

s2 = diff(smoothVecGauss(ampwav(s,winlen,winadv),1024/winadv));

% s2on = s2(1:round(length(s2)/3));
% pkInds = localMaxDer(abs(s2on),max(abs(s2on))/2);
% pkInds = pkInds(s2on(pkInds)>0);
% onset = pkInds(1);
% 
% s2off = s2(end-round(length(s2)/3):end);
% 
% pkInds = localMaxDer(abs(s2off),max(abs(s2off))/2);
% pkInds = pkInds(s2off(pkInds)<0);
% offset = pkInds(end) + round(length(s2)*2/3) - 1;

pkInds = localMaxDer(abs(s2),max(abs(s2))/3);
pkIndsOn = pkInds(s2(pkInds)>0);
onset = pkIndsOn(1);
pkIndsOff = pkInds(s2(pkInds)<0);
offset = pkIndsOff(end);


l = (offset - onset)*winadv*1000/44100;
