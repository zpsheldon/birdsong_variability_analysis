function h = day_division(wavdir,timesource)

if nargin==0
    wavdir = uigetdir(userpath,'pick directory');
end

if nargin<2
    timesource = questdlg('Get song times from file name or timestamps?','Song time source','File name','Timestamps','File name');
end

if exist([wavdir '\wavdirinfo.mat'],'file')
    load([wavdir '\wavdirinfo.mat'])
else
    dirstrct = mk_wavdirinfo(wavdir,timesource);
end

brkflg = 0;

for dayind = 1:length(dirstrct.dayset)
    
    dayind
    
    dayinds = find(dirstrct.daynums==dirstrct.dayset(dayind));
    dirnew = [wavdir '\day' num2str(dayind)];
    
    mkdir(dirnew);
    
    for flind = 1:length(dayinds)
        orgfile = [wavdir '\' dirstrct.wavfls{dayinds(flind)}];
        newfile = [dirnew '\' dirstrct.wavfls{dayinds(flind)}];
        flg = movefile(orgfile,dirnew);
        
        if ~flg
            errordlg('file failed to copy');
            brkflg = 1;
            break
        end
        
    end
    
    if brkflg
       break 
    end
    
    mk_wavdirinfo(dirnew,timesource);
    
end