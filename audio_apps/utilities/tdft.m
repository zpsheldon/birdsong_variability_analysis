function  [amps, phases, times, freqs, spec] = tdft(signal,winlen,winadv,fs)
%[amps, phases, times, freqs, spec] = tdft(signal,winlen,winadv,fs)

win = hanning(winlen);

signal_win = window_vc(signal,winlen,winadv);
signal_win = signal_win .* repmat(win,1,size(signal_win,2));

spec = fft(signal_win);

if rem(winlen,2),
    spec = spec(1:(winlen+1)/2,:);
else
    spec = spec(1:winlen/2+1,:);
end

amps = abs(spec);

if nargout > 1
    phases = angle(spec);
    times = 1000*(0:size(spec,2)-1)*winlen / fs;
    freqs = (0:size(spec,1)-1)*fs / winlen;
end