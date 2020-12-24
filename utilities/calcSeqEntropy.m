% function E = calcSeqEntropy(Xarr,Marr,Iarr)

clear
clc

load C:\Users\Chris\Documents\RESEARCH\SCHMIDTLAB\jon_stutters\OR_269\pre\dir\tstbd.mat
[transarr,Xarr,Marr,matchnms,Xarru,Marru,Iarr,Parr,Earr] = vc2N_transmat(vc,min_seq_nm,min_seq_len);

lvlnm = length(Xarr);

Earr2{1} = Earr{1};
Parr2{1} = Parr{1};

for i = 2:lvlnm
    Earr2{i} = [Earr2{i-1}(Iarr{i},:) Earr{i}];
    Parr2{i} = [Parr2{i-1}(Iarr{i},:) Parr{i}];
end

