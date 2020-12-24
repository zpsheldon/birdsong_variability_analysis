function [distvc,wavarr] = ampdistvc(wavdir,wavfls,clipons,clipsamps,filtcoeffs)

clipnm = length(clipons);

distvc = zeros(clipnm*(clipnm-1)/2,1);
k = 0;

for clipind = 1:clipnm    
    samp1 = clipons(clipind);
    samp2 = clipons(clipind) + clipsamps(clipind) - 1;
    flnm = [wavdir filesep wavfls{clipind}];
    
    wavtmp = filter(filtcoeffs(1,:),filtcoeffs(2,:),wavread(flnm,double([samp1,samp2])));
    amp1 = log(sqrt(mean(wavtmp.^2)));
    
    distvc(clipind) = amp1;
end