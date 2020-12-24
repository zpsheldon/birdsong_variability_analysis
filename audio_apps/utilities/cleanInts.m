function [cleaned,sngDelete,intVec] = cleanInts(ints,stdNum,opt)

if nargin<3
    opt = 0;
end
    
[sngVec,intVec] = findOutliers(ints,stdNum,opt);
sngDelete = unique(sngVec);
cleaned = ints(setdiff(1:size(ints,1),sngDelete),:);