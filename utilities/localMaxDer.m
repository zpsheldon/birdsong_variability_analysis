function [peakInds,peakWidths,peakHeights,peakRatios,onsetInds,offsetInds]  = localMaxDer(vec,thresh)

%[peakInds,widthsMid,peakRatios,peakHeights,onsetInds,offsetInds]  = localMaxDer(vec,thresh)

if min(vec) < 0
    vecAbs = abs(vec);
    negFlg = 1;
else
    negFlg = 0;
    vecAbs = vec;
end

if nargin == 1
    thresh = 0;
end

vecLen = length(vecAbs);

vecDiff = diff(vecAbs);
signDiff = sign(vecDiff);
peakVecRaw = diff(signDiff);
peakInds = find(peakVecRaw == -2) + 1;
peakNum = length(peakInds);

peakOffsetInds = find(peakVecRaw > 0) + 1;
peakOnsetInds = peakOffsetInds;

vecDiff2 = diff(vecDiff);
signDiff2 = sign(vecDiff2);
peakVecRaw2 = diff(signDiff2);
peakTopOffsetInds = find(peakVecRaw2 == 2) + 2;
peakTopOnsetInds = find(peakVecRaw2 == -2) + 2;

ratioCoeff = vecLen/(3*mean(vecAbs));
peakWidths = zeros(peakNum,1);
peakRatios = peakWidths;
peakHeights = peakWidths;
onsetInds = peakWidths;
offsetInds = peakWidths;

if length(peakInds) > 0

for i = 1:length(peakInds)

    pInd = peakInds(i);

    [onsetInds(i),offsetInds(i)] = findBounds(pInd,peakOnsetInds,peakOffsetInds,length(vecAbs));
    
    vecSpan = vecAbs(onsetInds(i):offsetInds(i));
    midLine = vecAbs(pInd) / 2;
    
    topOnsetCand = onsetInds(i):pInd;
    oDiff = vecAbs(topOnsetCand) - midLine;
    oInd = find(oDiff > 0);

    if isempty(oInd)
        onsetTopInds(i) = -1;
    else
        oIndInd = find(oDiff(oInd) == min(oDiff(oInd)));
        onsetTopInds(i) = topOnsetCand(oInd(oIndInd));
    end

    topOffsetCand = pInd:offsetInds(i);
    oDiff = vecAbs(topOffsetCand) - midLine;
    oInd = find(oDiff > 0);

    if isempty(oInd)
        offsetTopInds(i) = -1;
    else
        oIndInd = find(oDiff(oInd) == min(oDiff(oInd)));
        offsetTopInds(i) = topOffsetCand(oInd(oIndInd));
    end
    
    peakWidths(i) = offsetTopInds(i) - onsetTopInds(i) + 1;
    peakHeights(i) = vecAbs(pInd) - mean([vecAbs(offsetInds(i)),vecAbs(onsetInds(i))]);
    peakRatios(i) = ratioCoeff * peakHeights(i) / peakWidths(i);
    
end

if size(vecAbs,1) == 1
    vecAbs = vecAbs';
end

% elimInd = find(peakRatios < .5 | vecAbs(peakInds) < thresh | offsetTopInds' == -1 | onsetTopInds' == -1);
valInd = find(peakRatios > .5 & vecAbs(peakInds) > thresh & offsetTopInds'~=-1 & onsetTopInds' ~= -1);
peakInds = peakInds(valInd);
peakWidths = peakWidths(valInd);
peakRatios = peakRatios(valInd);
peakHeights = peakHeights(valInd);
onsetInds = onsetInds(valInd);
offsetInds = offsetInds(valInd);
% 
% 
% if ~isempty(elimInd)
%     for i = 1:length(elimInd)
%         
%         peakInds = deleteVec(peakInds,elimInd(i));
%         peakWidths = deleteVec(peakWidths,elimInd(i));
%         peakRatios = deleteVec(peakRatios,elimInd(i));
%         peakHeights = deleteVec(peakHeights,elimInd(i));
%         onsetInds = deleteVec(onsetInds,elimInd(i));
%         offsetInds = deleteVec(offsetInds,elimInd(i));
%         
%         elimInd = elimInd - 1;
%         
%     end
% end

onsetInds = onsetInds';
offsetInds = offsetInds';

end

%--------------------------------------------------------------------------
function [onset,offset] = findBounds(pInd,peakOnsetInds,peakOffsetInds,vecLen)

onsetIndDiff = peakOnsetInds - pInd;
indInd = find(onsetIndDiff < 0);

if ~isempty(indInd)
    onset = peakOnsetInds(max(indInd));
else
    onset = 1;
end

offsetIndDiff = peakOffsetInds - pInd;
indInd = find(offsetIndDiff > 0);

if ~isempty(indInd)
    offset = peakOffsetInds(min(indInd));
else
    offset = vecLen;
end
