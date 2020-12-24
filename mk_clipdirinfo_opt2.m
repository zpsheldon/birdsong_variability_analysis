function [clipstrct,dirstrct] = mk_clipdirinfo_opt2(wavdir)

if nargin==0
    wavdir = uigetdir(userpath,'pick directory');
end

if exist([wavdir filesep 'wavdirinfo.mat'])
    load([wavdir filesep 'wavdirinfo.mat'],'dirstrct')
else
    dirstrct = mk_wavdirinfo(wavdir);
end

wavnm = length(dirstrct.wavfls);

[dmy,fs] = wavread([wavdir filesep dirstrct.wavfls{1}]);

config_strct = ABconfig;

filt_f1 = config_strct.freqmin;
filt_f2 = config_strct.freqmax;
winlen = 2^round(log2(fs*config_strct.ampwin/1000));
winadv = winlen;

logthresh = config_strct.clip_amp_cutoff/winlen;
minsamps = uint32(floor(config_strct.clip_minlen * fs / 1000));
maxsamps = uint32(config_strct.clip_maxlen * fs / 1000);

[b,a] = butter(config_strct.filt_order,[filt_f1 filt_f2]/(fs/2));


gap_thresh_vc = 0:5:20;
gap_thresh_vc = gap_thresh_vc * fs / 1000;

paramnm = length(gap_thresh_vc);

zerostmp = zeros(1,10000,'uint32');

for paramind = 1:paramnm  
    clipinds{paramind} = zerostmp;
    wavinds{paramind} = zerostmp;
    clipsamps{paramind} = zerostmp;
    clipons{paramind} = zerostmp;  
    clipgaps{paramind} = zerostmp;  
end

h = waitbar(0/wavnm,['Making clip headers for ' wavdir]);

clipindvc = ones(1,paramnm);

for wavind = 1:wavnm
    
    if wavind>1
        h = waitbar(wavind/wavnm,h);
    end
    
    [s,fs] = wavread([wavdir filesep dirstrct.wavfls{wavind}]);
end