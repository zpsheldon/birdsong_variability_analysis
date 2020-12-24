function [gVal,gCrit] = calculateG(ints,varargin)

if ~isempty(varargin)
    biasInd = find(strcmp(varargin{:,1},'bias'));
    bias = varargin{biasInd,2};
else
    bias = 0;
end

seqIntNum = size(ints,2);

seqLens = sum(ints');
seqLenVar = var(seqLens);
intVar = sum(var(ints));

sampleNum = size(ints,1);
df = sampleNum - 1;
bias = bias / 2;

seqLenVar = seqLenVar - sampleNum*bias^2/df;

intVar = intVar - sampleNum*seqIntNum*bias^2/df;

gVal = seqLenVar/intVar;
%gVal = 1 - intVar/seqLenVar;

% if gVal < 0
%     gCrit = 1 - df/chi2inv(.025,df);
% else
%     gCrit = 1 - df/chi2inv(.975,df);
% end


% gVal = seqLenVar - intVar;
% 
if gVal < 0
    gCrit = intVar*(df/chi2inv(.975,df) - 1);
else
    gCrit = intVar*(df/chi2inv(.025,df) - 1);
end

