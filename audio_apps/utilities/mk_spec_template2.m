function [spmn,spstd,spalgn] = mk_spec_template2(specmat,sylstr)

sampsz = size(specmat,1);
lenvc = zeros(sampsz,1);

for i = 1:sampsz
    lenvc(i) = size(specarr{i},2);
end

freqnm = size(specarr{1},1);

specLen = median(lenvc);
mninds = find(lenvc==specLen);

if isempty(mninds)
    [dmy,mnind] = min(abs(lenvc-specLen));
    specLen = lenvc(mnind(1));
    mninds = find(lenvc==specLen);
end

if length(mninds)>1
    spalgn = zeros(length(mninds),freqnm,specLen);
    
    for i = 1:length(mninds)
        spalgn(i,:,:) = specarr{mninds(i)};
    end
    
    spmn = squeeze(mean(spalgn));
    
else
    spmn = specarr{mninds};
end

spalgn = zeros(sampsz,freqnm,specLen);

h = waitbar(1/sampsz,['Warping template ' sylstr]);

for sngInd = 1:sampsz
    
    spec = specarr{sngInd};
    w = dtwDist2(spmn,spec,[1.5 1.5 1]);
    
    [dmy,winds] = sort(w(:,1));
    w = w(winds,:);
    winds = find(w(:,1) > 0);
    
    spalgn(sngInd,:,w(winds,1)) = spec(:,round(w(winds,2)));
    
    h = waitbar(sngInd/sampsz,h,['Warping template ' sylstr]);
end

spmn = squeeze(mean(spalgn));
spstd = squeeze(std(spalgn));

close(h)