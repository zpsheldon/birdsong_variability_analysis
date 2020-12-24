function flg = check_headers(wavdir)

load([wavdir filesep 'wavdirinfo.mat'],'dirstrct')

fldir = dir(wavdir);
flnum = length(fldir)-2;
daynums = zeros(1,flnum);

flg = 0;

for flind = 1:flnum
    [dmy,flsnew{flind},ext{flind}] = fileparts(fldir(flind+2).name);
    flsnew{flind} = [flsnew{flind} '.wav'];
end

flsnew = flsnew(strcmp(ext,'.wav'));
flscommon = intersect(flsnew,dirstrct.wavfls);

if length(flscommon) ~= length(flsnew) || length(flscommon) ~= length(dirstrct.wavfls)
    flg = 1;
    inpt = questdlg('.wav files have been added or deleted since header file generation, regenerate?','Header warning','Yes','No','Yes');
    if strcmp(inpt,'Yes')
        flg = 2;
        mk_wavdirinfo(wavdir);
    end
    
    if flg==2 && exist([wavdir filesep 'clipdirinfo.mat'])
        inpt = questdlg('Regenerate clip header file?','Clip header warning','Yes','No','Yes');
        
        if strcmp(inpt,'Yes')
            mk_clipdirinfo(wavdir);
        end
    
    end
    
end

switch flg
    case 1
        htmp = warndlg('Deleted .wav files may cause errors, added files will be ignored');
    case 2
        htmp = warndlg('Header files have changed, it is recommended that any classification or analysis be restarted'); 
end

if exist('htmp')
    uiwait(htmp);
end