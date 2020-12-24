function [clipstrct,dirstrct] = mk_clipdirinfo_opt(wavdir)

if nargin==0
    wavdir = uigetdir(userpath,'pick directory');
end

if exist([wavdir filesep 'wavdirinfo.mat'])
    check_headers(wavdir);
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
    s = filter(b,a,s);
    
    for paramind = 1:paramnm
        
        gapthresh_tmp = gap_thresh_vc(paramind);
        
        [sampstmp,onstmp] = ...
            wav2clips(s,winlen,winadv,config_strct.amp_thresh,gapthresh_tmp,config_strct.normopt,config_strct.logopt,logthresh,minsamps,maxsamps);
        
        clipnm = length(onstmp);
        
        if clipnm > 1
            
            indstmp = clipindvc(paramind):clipindvc(paramind) + clipnm - 1;
            
            clipsamps{paramind}(indstmp) = sampstmp;
            clipons{paramind}(indstmp) = onstmp;
            
            offstmp = onstmp + sampstmp - 1;
            gapstmp = onstmp(2:end) - offstmp(1:end-1);
            
            wavinds{paramind}(indstmp) = wavind;
            clipinds{paramind}(indstmp) = 1:clipnm;
            
            clipindvc(paramind) = clipindvc(paramind) + clipnm; 
            
        end
        
    end
    
end

close(h)

clipindvc = clipindvc - 1;

for paramind = 1:paramnm  
    clipinds{paramind} = clipinds{paramind}(1:clipindvc(paramind));
    wavinds{paramind} = wavinds{paramind}(1:clipindvc(paramind));
    clipsamps{paramind} = clipsamps{paramind}(1:clipindvc(paramind));
    clipons{paramind} = clipons{paramind}(1:clipindvc(paramind));  
end

config_strct.fs = fs;
config_strct.b = b;
config_strct.a = a;

winind = 1;
pospropmat = zeros(paramnm-1,2);
Nmat = pospropmat;

h = waitbar(0/(paramnm-1),['Optimizing segmentation parameters for ' wavdir]);

for paramind = 2:paramnm
    
    ons1 = double(clipons{winind});
    wavinds1 = double(wavinds{winind});
    samps1 = double(clipsamps{winind});
    
    ons2 = double(clipons{paramind});
    wavinds2 = double(wavinds{paramind});
    samps2 = double(clipsamps{paramind});
    
    [posprop1,posprop2,N1,N2] = comparams(ons1,wavinds1,samps1,ons2,wavinds2,samps2,wavdir,config_strct,dirstrct);
    
    
    pospropmat(paramind-1,1) = posprop1;
    pospropmat(paramind-1,2) = posprop2;
    Nmat(paramind-1,1) = N1;
    Nmat(paramind-1,2) = N2;
    
    if posprop2 > posprop1
        winind = paramind;
    end
    
    h = waitbar((paramind-1)/(paramnm-1),h);
        
end


clipstrct.clipsamps = clipsamps{winind};
clipstrct.clipons = clipons{winind};
clipstrct.clipinds = clipinds{winind};
clipstrct.wavinds = wavinds{winind};
clipstrct.fs = fs;
clipstrct.gap_thresh = gap_thresh_vc(winind) * 1000 / fs;
clipstrct.amp_thresh = config_strct.amp_thresh;


clipstrct.filt_f1 = filt_f1;
clipstrct.filt_f2 = filt_f2;
clipstrct.filt_order = config_strct.filt_order;

clipstrct.gap_thresh_vc = gap_thresh_vc * 1000 / fs;
clipstrct.pospropmat = pospropmat;
clipstrct.Nmat = Nmat;

save([wavdir filesep 'clipdirinfo.mat'],'clipstrct')

close(h)

%--------------------------------------------------------------------------
function [posprop1,posprop2,N1,N2] = comparams(ons1,wavinds1,samps1,ons2,wavinds2,samps2,wavdir,config_strct,dirstrct)

[dmy,inds1,inds2] = setxor([ons1' wavinds1' samps1'],[ons2' wavinds2' samps2'],'rows');
N1 = length(inds1);
N2 = length(inds2);

sampnm = 50;
testnm = 1000;

sampinds1 = randperm(N1);
sampinds1 = inds1(sampinds1(1:(min(sampnm,length(sampinds1)))));

sampinds2 = randperm(N2);
sampinds2 = inds2(sampinds2(1:(min(sampnm,length(sampinds2)))));

testinds1 = setdiff(randperm(length(ons1)),sampinds1);
testinds1 = testinds1(1:min(length(testinds1),testnm));

testinds2 = setdiff(randperm(length(ons2)),sampinds2);
testinds2 = testinds2(1:min(length(testinds2),testnm));

distmat1 = zeros(length(testinds1),length(sampinds1));
distmat2 = zeros(length(testinds2),length(sampinds2));

freqvc = config_strct.fs*[0:round(512/2)]/512;
freqinds = find(freqvc >= config_strct.freqmin & freqvc <= config_strct.freqmax);

for sampind = 1:length(sampinds1)
    
    indtmp = sampinds1(sampind);
    samp1 = ons1(indtmp);
    samp2 = samp1 + samps1(indtmp) - 1;
    flnm = [wavdir filesep dirstrct.wavfls{wavinds1(indtmp)}];
    sampcliptmp = filter(config_strct.b,config_strct.a,wavread(flnm,double([samp1,samp2])));    
    
    sampcliptmp = tdft(sampcliptmp,512,512);
    sampclip{sampind} = log(max(sampcliptmp(freqinds,:),config_strct.amp_cutoff)) - log(config_strct.amp_cutoff);
    sampclip{sampind} = sampclip{sampind} / max(sampclip{sampind}(:));

end

for testind = 1:length(testinds1)
    indtmp = testinds1(testind);
    
    samp1 = ons1(indtmp);
    samp2 = samp1 + samps1(indtmp) - 1;
    flnm = [wavdir filesep dirstrct.wavfls{wavinds1(indtmp)}];
    testcliptmp = filter(config_strct.b,config_strct.a,wavread(flnm,double([samp1,samp2])));
    
    testcliptmp = tdft(testcliptmp,512,512);
    testclip{testind} = log(max(testcliptmp(freqinds,:),config_strct.amp_cutoff)) - log(config_strct.amp_cutoff);
    testclip{testind} = testclip{testind} / max(testclip{testind}(:));

end

for sampind = 1:length(sampinds1)
    for testind = 1:length(testinds1)        
        distmat1(testind,sampind) = dtwscore(sampclip{sampind},testclip{testind});
    end
end
    

sampclip = {};
testclip = {};

for sampind = 1:length(sampinds2)
    
    indtmp = sampinds2(sampind);
    samp1 = ons2(indtmp);
    samp2 = samp1 + samps2(indtmp) - 1;
    flnm = [wavdir filesep dirstrct.wavfls{wavinds2(indtmp)}];
    sampcliptmp = filter(config_strct.b,config_strct.a,wavread(flnm,double([samp1,samp2])));    
    
    sampcliptmp = tdft(sampcliptmp,512,512);
    sampclip{sampind} = log(max(sampcliptmp(freqinds,:),config_strct.amp_cutoff)) - log(config_strct.amp_cutoff);
    sampclip{sampind} = sampclip{sampind} / max(sampclip{sampind}(:));

end

for testind = 1:length(testinds2)
    indtmp = testinds2(testind);
    
    samp1 = ons2(indtmp);
    samp2 = samp1 + samps2(indtmp) - 1;
    flnm = [wavdir filesep dirstrct.wavfls{wavinds2(indtmp)}];
    testcliptmp = filter(config_strct.b,config_strct.a,wavread(flnm,double([samp1,samp2])));

    testcliptmp = tdft(testcliptmp,512,512);
    testclip{testind} = log(max(testcliptmp(freqinds,:),config_strct.amp_cutoff)) - log(config_strct.amp_cutoff);
    testclip{testind} = testclip{testind} / max(testclip{testind}(:));

end

for sampind = 1:length(sampinds2)
    for testind = 1:length(testinds2)        
        distmat2(testind,sampind) = dtwscore(sampclip{sampind},testclip{testind});
    end
end

posprop1 = length(find(distmat1(:) <= 1/25))/numel(distmat1);
posprop2 = length(find(distmat2(:) <= 1/25))/numel(distmat2);