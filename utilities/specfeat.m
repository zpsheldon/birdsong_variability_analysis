function featmat = specfeat(spec,freqs)

[fund,pwrat,ceps] = fund_freq(log(max(spec,.025)),freqs);
% w = wientropy(spec);
w = wientropy(ceps);
pwr = log(max(sqrt(sum(spec.^2)),.025));
featmat = [fund;w;pwr];