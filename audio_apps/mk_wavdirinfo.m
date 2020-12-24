function dirstrct = mk_wavdirinfo(wavdir,timesource)

if nargin==0
    wavdir = uigetdir(userpath,'pick directory');
end

fldir = dir(wavdir);
flnum = length(fldir)-2;
daynums = zeros(1,flnum);

if nargin<2
    timesource = questdlg('Get song times from file name or timestamps?','Song time source','File name','Timestamps','File name');
end

datevec = zeros(1,6);

for flind = 1:flnum
    [dmy,fls{flind},ext{flind}] = fileparts(fldir(flind+2).name);
    
    if strcmp(ext{flind},'.wav')
        if strcmp(timesource,'File name')
            datetmp = fls{flind}(end-18:end);
            dashinds = [0 findstr(datetmp,'_') length(datetmp)+1];
            
            for dashind = 1:6
                datevec(dashind) = str2num(datetmp(dashinds(dashind)+1:dashinds(dashind+1)-1));
            end
            
            daynums(flind) = datenum(datevec);
        else
            daynums(flind) = datenum(fldir(flind+2).date);
        end
        
        flsfull{flind} = [fls{flind} '.wav'];
    end
end

dirstrct.wavinds = find(strcmp(ext,'.wav'));
dirstrct.wavfls = flsfull(dirstrct.wavinds);
dirstrct.daynums = floor(daynums(dirstrct.wavinds));
dirstrct.dayset = unique(dirstrct.daynums);
dirstrct.daytms = daynums(dirstrct.wavinds);

save([wavdir filesep 'wavdirinfo.mat'],'dirstrct')

