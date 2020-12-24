clear
clc
close all

% WH 27

buff = 500;

load('/Users/cg/Documents/MATLAB/NE_data/WH27/timedata.mat')
load('/Users/cg/Documents/MATLAB/NE_data/WH27/UNOTHING/matchdata.mat')
load('/Users/cg/Documents/MATLAB/NE_data/WH27/UNOTHING/wavdirinfo.mat')

% load('/Users/cg/Documents/MATLAB/NE_data/WH27/FSAL_honed/matchdata.mat')
% load('/Users/cg/Documents/MATLAB/NE_data/WH27/FSAL_honed/wavdirinfo.mat')

strct = timestrct_U_sal;

seqN = numel(strct.sylarr);

[lenseq,sngDlt] = cleanInts(strct.lenseq,5,1);
Q = zeros(numel(strct.seqarr),0);
[W,psi,phi,sigma,omega] = CFAfull_prior_spc(lenseq,Q,200,1);

[z,u,eta] = timing_latvar(lenseq,W,sigma,psi);

indsval = setdiff(1:size(strct.lenseq,1),sngDlt);
wavinds = strct.wavinds(indsval);
seqinds = strct.seqinds(indsval);

lensUsal = sum((eta').^2)';

lensUsal = sum((eta').^1)';

lensUsal = z;

%lensUsal = sum(lenseqUsal')';
% lensFsal = sum(timestrct_F_sal.lenseq')';
% lensUNE = sum(timestrct_U_NE.lenseq')';

shortind = find(lensUsal==min(lensUsal));
longind = find(lensUsal==max(lensUsal));



shortwavind = wavinds(shortind);
longwavind = wavinds(longind);

shortmatchind = seqinds(shortind);
longmatchind = seqinds(longind);

[vc,labtmp,seqinds2,clipons2,cliplens2] = ...
    labs2vc(matchstrct.cliplabs,matchstrct.wavinds,matchstrct.cliptms,matchstrct.cliplens,...
    200);

seqonshort = clipons2(shortmatchind);
seqoffshort = clipons2(shortmatchind+seqN-1)+cliplens2(shortmatchind+seqN-1);

seqonlong = clipons2(longmatchind);
seqofflong = clipons2(longmatchind+seqN-1)+cliplens2(longmatchind+seqN-1);

wavflshort = dirstrct.wavfls{shortwavind};
wavfllong = dirstrct.wavfls{longwavind};

fs = 44100;

samponshort = seqonshort*fs/1000;
sampoffshort = seqoffshort*fs/1000;

samponlong = seqonlong*fs/1000;
sampofflong = seqofflong*fs/1000;


propstrct = ABconfig;
amp_cutoff = propstrct.amp_cutoff;

sshort = wavread([matchstrct.wavdir '/' wavflshort]);
slong = wavread([matchstrct.wavdir '/' wavfllong]);

sshort = sshort(samponshort:min(length(sshort),sampoffshort));
slong = slong(samponlong:min(length(slong),sampofflong));

[amps, phases, timeshort, freqs] = tdft(sshort,512,126,fs);
specshort = log(max(amps,amp_cutoff)) - log(amp_cutoff);

[amps, phases, timeslong, freqs] = tdft(slong,512,126,fs);
speclong= log(max(amps,amp_cutoff)) - log(amp_cutoff);

% specshort = specshort - max(specshort(:));
% speclong = speclong - max(speclong(:));

indsval = freqs>500 & freqs<9500;

figure
imagesc(timeslong,freqs(indsval),speclong(indsval,:))
hold on
imagesc(timeshort,freqs(indsval)+max(freqs(indsval))+20,specshort(indsval,:))
set(gca,'ydir','normal')
colormap jet

ylim([min(freqs(indsval)) 2*max(freqs(indsval))+20])

% seqstrct.matchinds = indsval;
% seqstrct.wavinds = seqinds2;
% seqstrct.cliptms = clipons2;
% seqstrct.cliplens = cliplens2;
% seqstrct.vc = vc;
% seqstrct.labsu = labtmp;
