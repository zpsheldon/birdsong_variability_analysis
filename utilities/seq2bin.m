function [X,sequ] = seq2bin(seq)

sequ = unique(seq);
d = length(sequ);
N = length(seq);

X = zeros(d,N);

if iscell(seq)
    for i = 1:d
        X(i,:) = strcmp(seq,sequ{i});
    end
else
    for i = 1:d
        X(i,:) = seq==sequ(i);
    end
end
