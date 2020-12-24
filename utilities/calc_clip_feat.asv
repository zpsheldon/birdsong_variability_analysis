function [cliplens,cliptms,clipwr,wavinds] = calc_clip_feat(maindir,dirstrct,clipstrct)

clipnm = length(clipstrct.clipons);
cliplens = zeros(clipnm,1);
clipwr = cliplens;
cliptms = cliplens;
wavinds = cliplens;

config_strct = ABconfig;

fs = 44100;

filt_f1 = config_strct.freqmin;
filt_f2 = config_strct.freqmax;
winlen = 2^round(log2(fs*config_strct.ampwin/1000));
winadv = winlen;
gap_thresh = config_strct.gap_thresh * fs / 1000;

logthresh = config_strct.clip_amp_cutoff/winlen;
minsamps = uint32(floor(config_strct.clip_minlen * fs / 1000));
maxsamps = uint32(config_strct.clip_maxlen * fs / 1000);

[b,a] = butter(config_strct.filt_order,[filt_f1 filt_f2]/(fs/2));

wavindsu = unique(clipstrct.wavinds);

clipind0 = 1;

for wavind = 1:length(wavindsu)
    
    flnm = [maindir filesep dirstrct.wavfls{wavindsu(wavind)}];
    signal = wavread(flnm);
    signal = filter(b,a,signal);

    [sampstmp,onstmp] = ...
        wav2clips(signal,winlen,winadv,config_strct.amp_thresh,gap_thresh,config_strct.normopt,config_strct.logopt,logthresh,minsamps,maxsamps,0);
    
    clipnmtmp = length(sampstmp);
    
    for clipind = 1:clipnmtmp
        samp1 = onstmp(clipind);
        samp2 = onstmp(clipind) + sampstmp(clipind) - 1;
        stmp = signal(samp1:samp2);
        clipwr(clipind0+clipind-1) = stmp'*stmp;
    end
    
    indstmp = clipind0:clipind0+clipnmtmp-1;
    
    clipwr(indstmp) = log(clipwr(indstmp) ./ double(sampstmp'));
    clipwr(indstmp) = clipwr(indstmp) - mean(clipwr(indstmp));
    
    cliplens(indstmp) = sampstmp;
    cliptms(indstmp) = onstmp;
    wavinds(indstmp) = wavindsu(wavind);
    
    clipind0 = clipind0 + clipnmtmp;
    
end

clipwr = clipwr(1:clipind0-1);
cliplens = cliplens(1:clipind0-1);
cliptms = cliptms(1:clipind0-1);

% clipwr = clipwr ./ cliplens;
cliplens = 1000 * cliplens / fs;
cliptms = 1000 * cliptms / fs;