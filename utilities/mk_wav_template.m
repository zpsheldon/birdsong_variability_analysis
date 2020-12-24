function [wavmn,wavstd,wmat] = mk_wav_template(wavmat,tp)

[sampsz,wavlen] = size(wavmat);

wavmn = mean(wavmat);
wavmnTmp = zeros(sampsz,wavlen);

pwt = [1.5 1.5 1];

h = waitbar(1/sampsz,['Making wav template']);

for wavind = 1:sampsz

    wav = wavmat(wavind,:);
    
    if tp==1
        w = dtwDist2(wavmn,wav,pwt);
    else
        w = dtwProd(wavmn,wav,pwt);
    end
    
    [dmy,winds] = sort(w(:,1));
    w = w(winds,:);
    winds = find(w(:,1) > 0);

    wavmnTmp(wavind,w(winds,1)) = wav(round(w(winds,2)));

    h = waitbar(wavind/sampsz,h,['Making wav template']);
end

wavmn = mean(wavmnTmp);
wavstd = std(wavmnTmp);

wmat = zeros(size(wavmat));

close(h)

h = waitbar(1/sampsz,['Warping wav']);

for wavind = 1:sampsz

    wav = wavmat(wavind,:);
    %     w = dtwDist2(wavmn,wav,[1.5 1.5 1]);
    if tp==1
        w = dtwDist2(wavmn,wav,pwt);
    else
        w = dtwProd(wavmn,wav,pwt);
    end
    
    [dmy,winds] = sort(w(:,1));
    wmat(wavind,w(winds,1)) = w(winds,2);
  
    h = waitbar(wavind/sampsz,h,['Warping wav']);
end

close(h)