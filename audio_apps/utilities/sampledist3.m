function [distvc,spmn,indOn,indOff,spalgn] = sampledist3(clipstrct,lab,N)
%CG: make return vectors within-cond and across-cond distances

indsval = find(strcmp(clipstrct.speclabs,lab));

if nargin==3
   indtmp = randperm(length(indsval));
   indsval = indsval(indtmp(1:N));
end

clipstrct.speclabs = clipstrct.speclabs(indsval);
clipstrct.wavarr = clipstrct.wavarr(indsval);

propstrct = ABconfig;

clipnm = length(clipstrct.speclabs);


propstrct.f_winlen = 2^round(log2(clipstrct.fs * propstrct.freqwin / 1000));

freqs = clipstrct.fs*[0:floor(propstrct.f_winlen/2)+1]/propstrct.f_winlen;
freqs = freqs(clipstrct.freqinds);

for clipind = 1:clipnm
    clip = clipstrct.wavarr{clipind};
    specarr{clipind} = tdft(clip,propstrct.f_winlen,propstrct.f_winlen/2);
    specarr{clipind} = log(max(specarr{clipind}(clipstrct.freqinds,:),propstrct.amp_cutoff));   
    specarr{clipind} = specarr{clipind} - log(propstrct.amp_cutoff);
end

[spmn,spstd,spalgn] = mk_spec_template(specarr,lab);

pkTms = getPks(spmn);
%     close(gcf)

indOn = min(pkTms);
indOff = max(pkTms);

distvc = zeros(clipnm,1);

h = waitbar(0/clipnm,'Calculating distances');

for clipind = 1:clipnm
    clip = specarr{clipind}; 
    distvc(clipind) = sqrt(mean((vec(spalgn(clipind,:,indOn:indOff)) - vec(spmn(:,indOn:indOff))).^2));
    h = waitbar(clipind/clipnm,h);
end

close(h)
