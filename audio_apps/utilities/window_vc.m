function s_win = window_vc(s,winlen,winadv)

s_len = length(s);
s = s(:); 
ncol = fix((s_len-winlen+winadv)/winadv);
% colind = 1 + (0:(ncol-1))*winadv;
% rowind = (1:winlen)';

s_win = s(repmat([1:winlen]',1,ncol) + repmat(winadv*(0:ncol-1),winlen,1));