function h = plotWavSpec(flnm)

if nargin == 0
    [filenm,pathnm] = uigetfile('*.wav','Pick a wav file');
    flnm = [pathnm filenm];
end

[s,fs] = wavread(flnm);

h = figure;

[amps, phases, times, freqs, spec] = tdft(s,512,256,fs);

logamps = log(max(amps,0.01));

imagesc(times,freqs,logamps);
set(gca,'ydir','normal');

ylim([300 9300])
