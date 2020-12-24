function featmat = specfeat2(spec,freqs)

[fund,pwrat,ceps] = fund_freq(spec,freqs);
% w = wientropy(spec);
w = wientropy(ceps);
pwr = sum(spec);
featmat = [fund;w;pwr];