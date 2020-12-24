function  s = itdft(spec,winlen,winadv)

[rows,cols] = size(spec);

winrat = winlen / winadv;
slen = (cols -1) * winadv + winlen;
s = zeros(winrat,slen);

if rem(winlen,2),
    spec = [spec;conj(flipud(spec(2:end,:)))];
else
    spec = [spec;conj(flipud(spec(2:end-1,:)))];
end

winmat = repmat(hanning(winlen),1,cols);

for i = 1:winrat
    stmp = real(ifft(spec(:,i:winrat:end))) .* winmat(:,i:winrat:end);
    sinds = (i-1)*winadv+1:(i-1)*winadv+length(stmp(:));
    s(i,sinds) = stmp(:)';
end

if winrat > 1
    s = mean(s);
end