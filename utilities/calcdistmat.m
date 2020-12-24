function dist = calcdistmat(specarr)

clipnm = length(specarr);

h = waitbar(0/clipnm,'Calculating distances');
dist = zeros(clipnm);

for clipind = 1:clipnm
    
    clip = specarr{clipind};
   
    for clipind2 = [1:clipind-1 clipind+1:clipnm]
        clip2 = specarr{clipind2};
        dist(clipind,clipind2) = dtwscore(clip,clip2);
    end
    
    h = waitbar(clipind/clipnm,h);
    
end

close(h)