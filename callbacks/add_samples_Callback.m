function add_samples_Callback(hObject,eventdata,handles)

clipstrct = get(handles.load_clips,'UserData');
propstrct = get(handles.clip_axis,'UserData');
load([clipstrct.wavdir filesep 'clipdirinfo.mat'])

clipstrct2 = clipstrct;

clipstrct = get(handles.load_clips,'UserData');
clipallnm = length(clipstrct2.clipsamps);
clipnm = length(clipstrct.speclabs);

sample_num = inputdlg('Sample number');
sample_num = str2num(sample_num{1});

nonm = clipallnm - clipnm;
nonsampinds = find(~ismember([clipstrct2.clipinds',clipstrct2.wavinds'],[clipstrct.clipinds',clipstrct.wavinds'],'rows'));

if ~isempty(sample_num)

    if sample_num > nonm
        sample_num = nonm;
        clipinds = nonsampinds;
        
        h = warndlg(['Sample number exceeds remaining clip number (' num2str(nonm) '), all selected'],'Warning');
        uiwait(h);

    else
        clipinds = randperm(nonm);
        clipinds = nonsampinds(clipinds(1:sample_num));
        
    end
    
else
    return;
end

clipnm = sample_num;

[dmy,indsrt] = sort(clipstrct2.clipsamps(clipinds));

clipstrct2.clipsamps = clipstrct2.clipsamps(clipinds(indsrt));
clipstrct2.wavinds = clipstrct2.wavinds(clipinds(indsrt));
clipstrct2.clipinds = clipstrct2.clipinds(clipinds(indsrt));
clipstrct2.clipons = clipstrct2.clipons(clipinds(indsrt));

freqs = propstrct.fs*[0:propstrct.f_winlen/2]/propstrct.f_winlen;

clipstrct2.speclens = zeros(1,clipnm);
[b,a] = butter(propstrct.filt_order,2*[propstrct.freqmin propstrct.freqmax]/propstrct.fs);

for clipind = 1:clipnm
    samp1 = clipstrct2.clipons(clipind);
    samp2 = clipstrct2.clipons(clipind) + clipstrct2.clipsamps(clipind) - 1;
    flnm = [clipstrct.wavdir filesep clipstrct.dirstrct.wavfls{clipstrct2.wavinds(clipind)}];
    clipstrct2.wavarr{clipind} = filter(b,a,wavread(flnm,double([samp1,samp2])));
    clipstrct2.specarr{clipind} = tdft(clipstrct2.wavarr{clipind},propstrct.f_winlen,propstrct.f_winadv);
    
    clipstrct2.featarr{clipind} = specfeat(clipstrct2.specarr{clipind}(clipstrct.freqinds,:),freqs(clipstrct.freqinds));
    
    clipstrct2.specarr{clipind} = log(max(clipstrct2.specarr{clipind}(clipstrct.freqinds,:),propstrct.amp_cutoff)) - log(propstrct.amp_cutoff);   
    clipstrct2.speclens(clipind) = size(clipstrct2.specarr{clipind},2);
    clipstrct2.speclabs{clipind} = ['x'];
end

clipstrct2.clipsamps = [clipstrct.clipsamps clipstrct2.clipsamps];
clipstrct2.wavinds = [clipstrct.wavinds clipstrct2.wavinds];
clipstrct2.clipinds = [clipstrct.clipinds clipstrct2.clipinds];
clipstrct2.clipons = [clipstrct.clipons clipstrct2.clipons];
clipstrct2.speclens = [clipstrct.speclens clipstrct2.speclens];
clipstrct2.wavarr = [clipstrct.wavarr clipstrct2.wavarr];
clipstrct2.specarr = [clipstrct.specarr clipstrct2.specarr];
clipstrct2.speclabs = [clipstrct.speclabs clipstrct2.speclabs];

clipnm = length(clipstrct2.wavinds);

clipstrct2.clipnm_tot = clipstrct.clipnm_tot;
clipstrct2.wavnm = clipstrct.wavnm;
clipstrct2.wavdir = clipstrct.wavdir;
clipstrct2.freqinds = clipstrct.freqinds;
clipstrct2.dirstrct = clipstrct.dirstrct;

clipstrct = clipstrct2;

propstrct.selectvc = zeros(1,clipnm);

distopt = questdlg('Compute clip distance matrix?','Distance prompt','Yes','No','Yes');
if strcmp(distopt,'Yes')
    for clipind = 1:clipnm
        clipstrct.specarr2{clipind} = tdft(clipstrct.wavarr{clipind},propstrct.f_winlen,propstrct.f_winlen);
        clipstrct.specarr2{clipind} = log(max(clipstrct.specarr2{clipind}(clipstrct.freqinds,:),propstrct.amp_cutoff)) - log(propstrct.amp_cutoff);
    end
    clipstrct.distmat = calcdistmat(clipstrct.specarr2);
else
    
    if isfield(clipstrct,'distmat')
        clipstrct = rmfield(clipstrct.distmat);
    end
    
end

set(handles.clip_axis,'UserData',propstrct);
set(handles.load_clips,'UserData',clipstrct);

set(handles.clip_slider,'Value',1)
plotspecstrct(handles,handles.clip_axis,handles.clip_slider,clipstrct,1000,5,'label');

set(handles.title_bar,'String',[clipstrct.wavdir ' (' num2str(clipstrct.wavnm) ' wav files, ' num2str(clipstrct.clipnm_tot) ' total clips, ' num2str(clipnm) ' samples)']);

