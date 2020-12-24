function [clipstrct,dirstrct] = mk_clipdirinfo(wavdir,timesource)

if nargin==0
    wavdir = uigetdir(userpath,'pick directory');
end

if exist([wavdir filesep 'wavdirinfo.mat'])
    load([wavdir filesep 'wavdirinfo.mat'],'dirstrct')
else
    if nargin<2
        dirstrct = mk_wavdirinfo(wavdir);
    else
        dirstrct = mk_wavdirinfo(wavdir,timesource);
    end
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
buff = 32*winadv;

[b,a] = butter(config_strct.filt_order,[filt_f1 filt_f2]/(fs/2));

clipinds = zeros(1,10000,'uint32');
wavinds = clipinds;
clipsamps = clipinds;
clipons = clipinds;

gap_thresh = config_strct.gap_thresh * fs / 1000;

clipind = 1;

h = waitbar(0/wavnm,'Making clip header');

for wavind = 1:wavnm
    
    if wavind>1
        h = waitbar(wavind/wavnm,h,'Making clip header');
    end
    
    [s,fs] = wavread([wavdir filesep dirstrct.wavfls{wavind}]);
    s = filter(b,a,s);
    [sampstmp,onstmp] = ...
        wav2clips(s,winlen,winadv,config_strct.amp_thresh,gap_thresh,config_strct.normopt,config_strct.logopt,logthresh,minsamps,maxsamps,buff);
    
    clipnm = length(onstmp);
    
    if clipnm > 1
        
        indstmp = clipind:clipind + clipnm - 1;
        
        clipsamps(indstmp) = sampstmp;
        clipons(indstmp) = onstmp;
        wavinds(indstmp) = wavind;
        clipinds(indstmp) = 1:clipnm;
        
        clipind = clipind + clipnm;
        
    end
    
end


clipnm = clipind - 1;

clipstrct.clipsamps = clipsamps(1:clipnm);
clipstrct.clipons = clipons(1:clipnm);
clipstrct.clipinds = clipinds(1:clipnm);
clipstrct.wavinds = wavinds(1:clipnm);
clipstrct.fs = fs;


clipstrct.filt_f1 = filt_f1;
clipstrct.filt_f2 = filt_f2;
clipstrct.filt_order = config_strct.filt_order;

save([wavdir filesep 'clipdirinfo.mat'],'clipstrct')

close(h)