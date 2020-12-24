function [outMat,colStrct] = findOutliers(data,stdNum,opt)

[rowNum,colNum] = size(data);

if opt == 0
    threshVec = std(data)*stdNum;
    dataDev = abs(data - ones(rowNum,1)*mean(data));
else
    threshVec = iqr(data)*stdNum*0.7413;
    dataDev = abs(data - ones(rowNum,1)*median(data));
end

devSum = sum(dataDev');
best = find(devSum == min(devSum));

outMat = [];

for i = 1:colNum
    thresh = threshVec(i);
    inds = find(dataDev(:,i) >= thresh);
    outMat = [outMat; inds];
    colStrct{i} = inds;
end