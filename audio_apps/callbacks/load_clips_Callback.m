function load_clips_Callback(hObject,eventdata,handles)

clipdir = uigetdir(userpath,'pick directory');

if ~clipdir
    return
end

if ~(exist([clipdir filesep 'clipdirinfo.mat']))
    button = questdlg('No clip header file, make it?');
    if strcmp(button,'Yes')
        [clipstrct,dirstrct] = mk_clipdirinfo(clipdir);
    else
        return;
    end
else
    
    flg = check_headers(clipdir);
    
    load([clipdir filesep 'clipdirinfo.mat'])
    load([clipdir filesep 'wavdirinfo.mat'])
end

set(handles.title_bar,'String',clipdir);

dirstrct.wavdir = clipdir;
clipstrct.wavdir = clipdir;

clipstrct.dirstrct = dirstrct;

propstrct = get(handles.clip_axis,'UserData');
propstrct.f_winlen = 2^round(log2(clipstrct.fs * propstrct.freqwin / 1000));
propstrct.f_winadv = 2^round(log2(clipstrct.fs * propstrct.timewin / 1000));
propstrct.fs = clipstrct.fs;

clipnm = length(clipstrct.clipsamps);

sample_num = inputdlg('Sample number');
sample_num = str2num(sample_num{1});

if ~isempty(sample_num)
    
    if sample_num > clipnm
        sample_num = clipnm;
        clipinds = 1:clipnm;
        
        h = warndlg(['Sample number exceeds clip number (' num2str(clipnm) '), all selected'],'Warning');
        uiwait(h);
        
    else
        clipinds = randperm(clipnm);
        clipinds = clipinds(1:sample_num);
        
    end
    
    
else
    sample_num = clipnm;
    clipinds = 1:clipnm;
end

clipnm_tot = clipnm;
clipnm = sample_num;

[dmy,indsrt] = sort(clipstrct.clipsamps(clipinds));

clipstrct.clipsamps = clipstrct.clipsamps(clipinds(indsrt));
clipstrct.wavinds = clipstrct.wavinds(clipinds(indsrt));
clipstrct.clipinds = clipstrct.clipinds(clipinds(indsrt));
clipstrct.clipons = clipstrct.clipons(clipinds(indsrt));

propstrct.selectvc = zeros(1,clipnm);

freqs = propstrct.fs*[0:propstrct.f_winlen/2]/propstrct.f_winlen;
clipstrct.freqinds = find(freqs >= propstrct.freqmin & freqs <= propstrct.freqmax);

clipstrct.speclens = zeros(1,clipnm);
[b,a] = butter(propstrct.filt_order,2*[propstrct.freqmin propstrct.freqmax]/propstrct.fs);

wavset = unique(clipstrct.wavinds);

for wavind = 1:length(wavset)
    
    clipinds = find(clipstrct.wavinds==wavset(wavind));
    flnm = [clipstrct.wavdir filesep clipstrct.dirstrct.wavfls{wavset(wavind)}];
    
    if exist(flnm)
        signal = wavread(flnm);
        
        for clipindsub = 1:length(clipinds)
            clipind = clipinds(clipindsub);
            samp1 = clipstrct.clipons(clipind);
            samp2 = clipstrct.clipons(clipind) + clipstrct.clipsamps(clipind) - 1;
            
            clipstrct.wavarr{clipind} = filter(b,a,signal(samp1:samp2));
            clipstrct.specarr{clipind} = tdft(clipstrct.wavarr{clipind},propstrct.f_winlen,propstrct.f_winadv);
            
            clipstrct.featarr{clipind} = specfeat(clipstrct.specarr{clipind},freqs(clipstrct.freqinds));
            
            clipstrct.specarr{clipind} = log(max(clipstrct.specarr{clipind}(clipstrct.freqinds,:),propstrct.amp_cutoff)) - log(propstrct.amp_cutoff);
            clipstrct.speclens(clipind) = size(clipstrct.specarr{clipind},2);
            clipstrct.speclabs{clipind} = ['x'];
        end
        
    else
        errordlg([flnm ' missing, sampling canceled'])
        return;
    end
end

wavnm = length(dirstrct.wavinds);
clipstrct.clipnm_tot = clipnm_tot;
clipstrct.wavnm = wavnm;

distopt = questdlg('Compute clip distance matrix?','Distance prompt','Yes','No','Yes');
if strcmp(distopt,'Yes')
    for clipind = 1:clipnm
        clipstrct.specarr2{clipind} = tdft(clipstrct.wavarr{clipind},propstrct.f_winlen,propstrct.f_winlen);
        clipstrct.specarr2{clipind} = log(max(clipstrct.specarr2{clipind}(clipstrct.freqinds,:),propstrct.amp_cutoff)) - log(propstrct.amp_cutoff);
    end
    clipstrct.distmat = calcdistmat(clipstrct.specarr2);
end

set(handles.clip_axis,'UserData',propstrct);
set(handles.load_clips,'UserData',clipstrct);

set(handles.clip_slider,'Value',1)
plotspecstrct(handles,handles.clip_axis,handles.clip_slider,clipstrct,1000,5,'label');

set(handles.title_bar,'String',[clipdir ' (' num2str(wavnm) ' wav files, ' num2str(clipnm_tot) ' total clips, ' num2str(clipnm) ' samples)']);
