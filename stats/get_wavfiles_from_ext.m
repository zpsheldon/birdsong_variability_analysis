clear
clc
close all

dirstrct = dir(['/Volumes/MELS SONG D/Steve/Infusion/WH27/UNOTHING/']);
savedir = ['WH27/UNOTHING/'];

flnum = numel(dirstrct)-2;


for flind = 1:1200 %flnum
    
    flind
    
    flnm = dirstrct(flind+2).name;
    gotit = 0;
    
    try
        [s,fs] = audioread(['/Volumes/MELS SONG D/Steve/Infusion/WH27/UNOTHING/' flnm]);
        
        if ~isempty(s)
            gotit = 1;
        end
        
    catch
        
    end
    
    if ~gotit
        flnmnew = ['a' flnm];
        % wavwrite(s,fs,[savedir flnm]);
        movefile(['/Volumes/MELS SONG D/Steve/Infusion/WH27/UNOTHING/' flnm],['/Volumes/MELS SONG D/Steve/Infusion/WH27/UNOTHING/' flnmnew]);
    end
    
end