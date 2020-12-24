function [P_seq,E_seq,cumN_seq,P_seqb,E_seqb,cumN_seqb] = calctranstrct(vc,seq,Nmin)

% P_seq = zeros(1,length(seq)-1);
% E_seq = P_seq;
% cumN_seq = P_seq;
% P_seqb = P_seq;
% E_seqb = P_seq;
% cumN_seqb = P_seq;

seqlen = length(seq);

[Xarr,Narr,Iarr,Parr,Earr] = vc2N_transmat_2seq(vc,seq,0,Nmin,0);

P_seq = cell2mat(Parr);
E_seq = cell2mat(Earr);
cumN_seq = cell2mat(Narr);

[Xarr,Narr,Iarr,Parr,Earr] = vc2N_transmat_2seq(vc,seq,0,Nmin,1);

P_seqb = fliplr(cell2mat(Parr));
E_seqb = fliplr(cell2mat(Earr));
cumN_seqb = fliplr(cell2mat(Narr));