function match_all_clips_Callback(hObject, eventdata, handles)
% hObject    handle to match_all_clips (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


templatestrct = get(handles.make_template,'UserData');
propstrct = get(handles.template_axis,'UserData');
trainstrct = get(handles.load_clips,'UserData');

usecur = 0;
if ~isempty(trainstrct)
    wavdir = trainstrct.wavdir;
    opt = questdlg('Use current clip data?','Data prompt','Yes','No','Yes');
    if strcmp(opt,'Yes')
        usecur = 1;
    end
end

if ~usecur
    [flnm,flpth] = uigetfile('*.mat','Choose training data');
    training_filename = [flpth flnm];
    if ~training_filename
        return;
    end
    
    load(training_filename);
    trainstrct = clipstrct;
    trainstrct.clipfl = training_filename;
    
    if ~exist('wavdir')
        wavdir = trainstrct.wavdir;
    end
    
end

flg = check_headers(wavdir);

load([wavdir filesep 'wavdirinfo.mat']);
load([wavdir filesep 'clipdirinfo.mat']);

freqinds = trainstrct.freqinds;
labsu = templatestrct.speclabs;
labsu{end+1} = 'x';
clipnm = length(clipstrct.clipons);

[flnm,flpth] = uiputfile('*.mat','Choose where to save match data');
match_filename = [flpth flnm];
if ~match_filename
    return;
end


clipnm = length(clipstrct.clipons);
useans = questdlg('Use manually labeled training data or match everything?','Training option','Use labeled','Match everything','Use labeled');
class = cell(1,clipnm);

switch useans
    case 'Use labeled'
        
        protonm = length(trainstrct.specarr);
        
        traininds = zeros(1,protonm);
        for protoind = 1:protonm
            traininds(protoind) = find(clipstrct.wavinds==trainstrct.wavinds(protoind) & clipstrct.clipinds==trainstrct.clipinds(protoind));
        end
        
        novinds = setdiff(1:clipnm,traininds);
        novnm = length(novinds);
        
        class(traininds) = trainstrct.speclabs;
        
        
    case 'Match everything'
        novinds = 1:clipnm;
        novnm = clipnm;
    otherwise
        return;
end


[b,a] = butter(propstrct.filt_order,2*[propstrct.freqmin propstrct.freqmax]/propstrct.fs);
winlen = propstrct.f_winlen;
winadv = propstrct.f_winlen;


f = propstrct.fs*[0:winlen/2]/winlen;
f = f(trainstrct.freqinds);

strct = ABconfig;
amp_cutoff = propstrct.amp_cutoff;

protonm = length(trainstrct.specarr);

clipstr = num2str(novnm);
msg = ['Matching clips (' clipstr ')'];

protodist = zeros(1,protonm);

h = waitbar(1/novnm,msg);

maxlen = ceil(max(clipstrct.clipsamps) / winadv);

zeromat = zeros(length(f),maxlen);

labnm = length(labsu);
clipcnt = zeros(labnm,1);

for clipind = 1:novnm
    
    samp1 = clipstrct.clipons(novinds(clipind));
    samp2 = clipstrct.clipons(novinds(clipind)) + clipstrct.clipsamps(novinds(clipind)) - 1;
    flnm = [wavdir filesep dirstrct.wavfls{clipstrct.wavinds(novinds(clipind))}];
    
    if exist(flnm)
        s = filter(b,a,wavread(flnm,double([samp1,samp2])));
        clip = tdft(s,winlen,winadv);
        clear s
        clip = log(max(clip(freqinds,:),amp_cutoff)) - log(amp_cutoff);
        
        for protoind = 1:protonm
            
            d = size(clip,2) - size(trainstrct.specarr2{protoind},2);
            
            if d > 0
                spec2 = [zeromat(:,1:abs(round(d/2))) trainstrct.specarr2{protoind} zeromat(:,1:abs(d)-abs(round(d/2)))];
                spec1 = clip;
            elseif d < 0
                spec1 = [zeromat(:,1:abs(round(d/2))) clip zeromat(:,1:abs(d)-abs(round(d/2)))];
                spec2 =  trainstrct.specarr2{protoind};
            else
                spec1 = clip;
                spec2 = trainstrct.specarr2{protoind};
            end
            
            protodist(protoind) = min(xcorr21(spec1,spec2,1));
            
        end
        
        [dmy,srtinds] = sort(protodist,2,'ascend');
        srtinds = srtinds(1:trainstrct.kopt);
        
        for labind = 1:labnm
           clipcnt(labind) = length(find(strcmp(trainstrct.speclabs(srtinds),labsu{labind})));
        end
        
        [maxnm,winind] = max(clipcnt);
        
        if length(winind)>1
            distmp = protodist(srtinds);
            winindind = 1;
            winscr = mean(distmp(find(strcmp(labsu{winind(winindind)},trainstrct.speclabs(srtinds)))));
            
            for indtmp = 2:length(winind)
                scrtmp = mean(distmp(find(strcmp(labsu{winind(winindind)},trainstrct.speclabs(srtinds)))));
                if scrtmp < winscr
                    winindind = indtmp;
                    winscr = scrtmp;
                end
            end
            
            winind = winind(winindind);
        end
        
        votes = clipcnt(winind);
       
        if ~strcmp(labsu{winind},'x') && votes >= trainstrct.threshvc(winind)
            class{novinds(clipind)} = labsu{winind};
        else
            class{novinds(clipind)} = 'x';
        end
        
        if rem(clipind,100)==0
            h = waitbar(clipind/novnm,h);
        end
        
    else
        errordlg([flnm ' missing, matching canceled'])
        return;
        
    end
    
    
end



close(h)

cliplens = 1000 * double(clipstrct.clipsamps) / propstrct.fs;
cliptms = 1000 * double(clipstrct.clipons) / propstrct.fs;

matchstrct.cliplabs = class;
matchstrct.wavdir = wavdir;
matchstrct.wavinds = clipstrct.wavinds;
matchstrct.clipinds = clipstrct.clipinds;
matchstrct.cliptms = cliptms;
matchstrct.cliplens = cliplens;
matchstrct.trainingfl = trainstrct.clipfl;

save(match_filename,'matchstrct','novinds')
