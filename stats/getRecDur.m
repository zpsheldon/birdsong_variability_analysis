function [recdur,wavindsval] = getRecDur(wavtms,daylims,maxgap,cutoff)

if nargin<4
    cutoff = inf;
end

recdur = 0;

wavtms = wavtms(:);

wavtms = wavtms(wavtms>daylims(1) & wavtms<daylims(2));

[wavtms,wavindsort] = sort(wavtms,'ascend');

t0 = wavtms(1);
begint0 = t0;

wavindsval = zeros(numel(wavtms),1);
wavindsval(1) = 1;

for wavind = 1:numel(wavtms)-1
    tnext = wavtms(wavind+1);
    
    if tnext<(t0+maxgap) & (tnext-begint0)<cutoff
        recdur = recdur + tnext - t0;
    end
    
    if tnext>=(t0+maxgap)
        begint0 = tnext;
    end
    
    if (tnext-begint0)<cutoff
       wavindsval(wavind+1)=1; 
    end
    
    t0 = tnext;
    
end

wavindsval = wavindsort(logical(wavindsval));