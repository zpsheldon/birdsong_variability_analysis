function distvc = specdistvc(wavdir,wavfls,clipons,clipsamps,winlen,winadv,filtcoeffs,freqinds,amp_cutoff)

clipnm = length(clipons);

distvc = zeros(clipnm*(clipnm-1)/2,1);
k = 0;

for clipind = 1:clipnm    
    samp1 = clipons(clipind);
    samp2 = clipons(clipind) + clipsamps(clipind) - 1;
    flnm = [wavdir filesep wavfls{clipind}];
    
    wavtmp = filter(filtcoeffs(1,:),filtcoeffs(2,:),wavread(flnm,double([samp1,samp2])));
    spec1 = tdft(wavtmp,winlen,winadv);
    
    spec1 = spec1 / median(spec1(:));
    
%     logspec1 = log(max(spec1(freqinds,:),amp_cutoff));
%     spec1 = spec1 - mean(spec1(:));
    
    for clipind2 = clipind+1:clipnm
        
        samp1 = clipons(clipind2);
        samp2 = clipons(clipind2) + clipsamps(clipind2) - 1;
        flnm = [wavdir filesep wavfls{clipind2}];
        
        wavtmp = filter(filtcoeffs(1,:),filtcoeffs(2,:),wavread(flnm,double([samp1,samp2])));
        spec2 = tdft(wavtmp,winlen,winadv);
        
        spec2 = spec2 / median(spec2(:));
        
%         logspec2 = log(max(spec2(freqinds,:),amp_cutoff));
        
%         spec2 = spec2 - mean(spec2(:));
        
        [w,D,d,scr]=dtwDist(spec1,spec2,[1.5 1.5 1]);
%         spec2hat = spec2(:,max(round(sort(w(:,2))),1));
%         spec1hat = spec1(:,max(round(sort(w(:,1))),1));
%         
% %         inds = find(spec1hat>log(amp_cutoff) & spec2hat > log(amp_cutoff));
% %         
%         b = robustfit(spec1hat(:),spec2hat(:),'bisquare',4.685,'off');
%         r = spec1hat - spec2hat/b;
%         
%         scr = sqrt(mean(r(:).^2));
% 
%         
%         
        k = k + 1;
 
        distvc(k) = scr;
        
    end
    
end