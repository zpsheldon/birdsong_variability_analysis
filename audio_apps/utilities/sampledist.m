function [classmat,classarr,distwithin,distbw,distmat] = sampledist(clipstrct1,clipstrct2,lab1,lab2)
%CG: make return vectors within-cond and across-cond distances

if nargin == 0
    tp = 'spectral';
end

if nargin < 2
    [filename1, pathname1] = uigetfile('*.mat', 'Pick 1st sample file');
    load([pathname1 filename1],'clipstrct')
    clipstrct1 = clipstrct;
    clear clipstrct
end


if nargin < 3
    [filename2, pathname2] = uigetfile('*.mat', 'Pick 2nd sample file');
    load([pathname2 filename2],'clipstrct')
    clipstrct2 = clipstrct;
    clear clipstrct
end


if nargin < 4
    lab1 = inputdlg('1st class');
    lab1 = lab1{1};
    lab2 = inputdlg('2nd class');
    lab2 = lab2{1};
end


indsval1 = strcmp(clipstrct1.speclabs,lab1);
indsval2 = strcmp(clipstrct2.speclabs,lab2);

clipstrct1.speclabs = clipstrct1.speclabs(indsval1);
clipstrct1.wavarr = clipstrct1.wavarr(indsval1);
clipstrct2.speclabs = clipstrct2.speclabs(indsval2);
clipstrct2.wavarr = clipstrct2.wavarr(indsval2);

propstrct = ABconfig;

clipnm1 = length(clipstrct1.speclabs);
clipnm2 = length(clipstrct2.speclabs);

clipnm = clipnm1 + clipnm2;
sampvc = [ones(clipnm1,1);2*ones(clipnm2,1)];

distmat = zeros(clipnm);

labsu1 = unique(clipstrct1.speclabs);
labsu2 = unique(clipstrct2.speclabs);

labstot = length(labsu1) + length(labsu2);

labvc = zeros(clipnm,1);


propstrct.f_winlen = 2^round(log2(clipstrct1.fs * propstrct.freqwin / 1000));

freqs = clipstrct1.fs*[0:floor(propstrct.f_winlen/2)+1]/propstrct.f_winlen;
freqs = freqs(clipstrct1.freqinds);

featmat = zeros(clipnm,3);

for labind = 1:length(labsu1)
    labsu{labind} = ['1-' labsu1{labind}];
end

for labind = 1:length(labsu2)
    labsu{labind+length(labsu1)} = ['2-' labsu2{labind}];
end


for clipind = 1:clipnm    
    if sampvc(clipind)==1
        indtmp = find(strcmp(labsu1,clipstrct1.speclabs{clipind}));
        labtmp = clipstrct1.speclabs{clipind};
        clip = clipstrct1.wavarr{clipind};
    else
        indtmp = find(strcmp(labsu2,clipstrct2.speclabs{clipind-clipnm1}))+length(labsu1);
        labtmp = clipstrct2.speclabs{clipind-clipnm1};
        clip = clipstrct2.wavarr{clipind-clipnm1};
    end
    
    labvc(clipind) = indtmp;
    labsall{clipind} = [num2str(sampvc(clipind)) '-' labtmp];
    specarr{clipind} = tdft(clip,propstrct.f_winlen,propstrct.f_winlen);
    specarr{clipind} = log(max(specarr{clipind}(clipstrct1.freqinds,:),propstrct.amp_cutoff));
    
end




matchtot = clipnm * (clipnm - 1) / 2;
k = 0;

h = waitbar(0/matchtot,'Calculating distance matrix');

for clipind = 1:clipnm
    
    clip1 = specarr{clipind};
    
    for clipind2 = clipind + 1:clipnm
        
        clip2 = specarr{clipind2};
        
        distmat(clipind,clipind2) = dtwscore2(clip1,clip2);
        distmat(clipind2,clipind) = distmat(clipind,clipind2);
        
        k = k + 1;
        
        h = waitbar(k/matchtot,h);
        
    end
    
    
end

close(h)

classmat = zeros(labstot);
for labind = 1:length(labsu1)
    labsu{labind} = ['1-' labsu1{labind}];
end

for labind = 1:length(labsu2)
    labsu{labind+length(labsu1)} = ['2-' labsu2{labind}];
end

indmat = repmat(1:clipnm,clipnm,1);
labmat = repmat(labvc,1,clipnm);

for labind = 1:labstot
    for labind2 = labind:labstot
        inds = find(indmat < indmat' & ((labmat == labind & labmat' == labind2) | (labmat' == labind & labmat == labind2)));
        classmat(labind,labind2) = median(distmat(inds));
        classmat(labind2,labind) = classmat(labind,labind2);
        classarr{labind,labind2} = distmat(inds);
    end
end

distwithin = [classarr{1,1};classarr{2,2}];
distbw = classarr{1,2};