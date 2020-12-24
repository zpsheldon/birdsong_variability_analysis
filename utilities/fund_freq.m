function [fund,pwrat,ceps,cep_freqs] = fund_freq(spec,f) 

[freqnm,tnm] = size(spec);

fund = zeros(1,tnm);
freqnm2 = 16*129;
% freqnm2 = 4*129;
freqmax = max(f)/2;
freqmax = 3000;
freqmin = 300;

fband = f(2)-f(1);
fsband = 1/fband;
freqs2 = fsband*[0:ceil(freqnm2/2)-1]/freqnm2;
freqs3 = 1./freqs2;

spec = spec - repmat(mean(spec),freqnm,1);
  
ceps = abs(fft(spec,freqnm2));
ceps = ceps(1:length(freqs2),:);

freqeval = find(freqs3 <= freqmax & freqs3 >= freqmin);
ceps = ceps(freqeval,:);
freqs3 = freqs3(freqeval);

[pwr,maxind] = max(ceps);
fund = freqs3(maxind);
pwrat = pwr ./ sum(ceps);

cep_freqs = freqs3;

% for tind = 1:tnm
% 
% 
%     [ceps,Pxxc,f] = pmtm(spec(:,tind),1.5,freqnm2);
%     maxind = max(2,find(ceps==max(ceps)));
%     
%     fund(tind) = 1/freqs2(maxind);
% 
%     
% end
