function [stdvc,distmat,stdtot,specmat,specdevmatnorm,timeinds] = spec_regression(spalgnarr)

sylN = numel(spalgnarr);

spectmp = squeeze(spalgnarr{1}(1,:,:));
freqN = size(spectmp,1);

timeN = 1;

timeinds = cell(sylN,1);

for sylind = 1:sylN
    timeinds{sylind} = timeN:timeN+size(spalgnarr{sylind},3)-1;
    timeN = timeN + size(spalgnarr{sylind},3);
end

seqN = size(spalgnarr{1},1);

specdevmat = zeros(seqN,freqN,timeN);
specdevmatnorm = specdevmat;
specmat = specdevmat;

for sylind = 1:sylN
    spmn = squeeze(mean(spalgnarr{sylind}));
    for seqind = 1:seqN
        specmat(seqind,:,timeinds{sylind})= squeeze(spalgnarr{sylind}(seqind,:,:));
        specdevmat(seqind,:,timeinds{sylind})= squeeze(spalgnarr{sylind}(seqind,:,:))-spmn;
    end
end

for seqind = 1:seqN
    freqdev = mean(specdevmat(seqind,:,:),3)';
    freqdevmat = repmat(freqdev,1,timeN);
    specdevmatnorm(seqind,:,:)= squeeze(specdevmat(seqind,:,:)) - freqdevmat;
end

distmat = zeros(seqN,sylN);

for sylind = 1:sylN
    for seqind = 1:seqN
        distmat(seqind,sylind) = mean(vec(specdevmatnorm(seqind,:,timeinds{sylind}).^2));
    end
end

stdvc = sqrt(mean(distmat));
stdtot = std(specdevmatnorm(:));


