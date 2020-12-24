function amps = ampwav(signal,winlen,winadv)

signal_win = window_vc(signal,winlen,winadv);
amps = sqrt(sum(signal_win .^ 2) / winlen);