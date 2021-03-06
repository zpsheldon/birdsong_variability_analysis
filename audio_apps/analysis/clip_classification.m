function varargout = clip_classification(varargin)

h_mainfig = clipclassification_mainfig();

handles = guihandles(h_mainfig);

handlenms = fieldnames(handles);

handles.clickfun = @clipick;

for i = 1:length(handlenms)
    h = getfield(handles,handlenms{i});
    if isfield(get(h),'Callback') & strcmp(get(h,'Callback'),'1')
        eval(['set(h,''Callback'',','{@' get(h,'Tag') '_Callback,handles});'])
    end
        
end

% set(handles.load_clips,'Callback',{@load_clips_Callback,handles});
set(handles.clip_slider,'Value',1);

propstrct = ABconfig;

set(handles.clip_axis,'UserData',propstrct);
set(handles.template_axis,'UserData',propstrct);

hcolor(handles);

h = subplot(1,1,1,'Parent',handles.analysis_panel);
set(h,'ytick',[],'xtick',[])

% --- Executes on button press in make_template.
function make_template_Callback(hObject, eventdata, handles)
% hObject    handle to make_template (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


templatestrct = get(handles.make_template,'UserData');
clip_propstrct = get(handles.clip_axis,'UserData');
clipstrct = get(handles.load_clips,'UserData');
propstrct = get(handles.template_axis,'UserData');


if ~isfield(templatestrct,'freqinds')
    templatestrct.freqinds = clipstrct.freqinds;
end

if ~isfield(templatestrct,'wavdir')
    templatestrct.wavdir = clipstrct.wavdir;
end

classnm = inputdlg('Class label');

if strcmp(classnm,'x')
    errordlg('''x'' reserved for non-matches, choose another label')
    return;
end

template_ind = 1;
exemplar_nm = 0;

% template structure consists of 'speclabs', 'template_inds',
% 'template'

if isfield(templatestrct,'speclabs')
    template_ind = find(strcmp(templatestrct.speclabs,classnm));
else
    templatestrct.speclabs = {};
    templatestrct.threshvc = [];
    template_ind = [];
end

clipstmp = clipstrct.specarr(find(clip_propstrct.selectvc));
clipnm = length(clipstmp);

if isempty(template_ind)
    template_ind = length(templatestrct.speclabs)+1;
end

templatestrct.speclabs(template_ind) = classnm;
templatestrct.threshvc(template_ind) = 20;

% templatestrct.specarr{template_ind} = mkTemplate(clipstrct.specarr(templatestrct.specinds{template_ind}),0);
clip_selectinds = find(clip_propstrct.selectvc);
if length(clip_selectinds)>1
    templatestrct.specarr{template_ind} = mkTemplateDTW(clipstrct.specarr(clip_selectinds));
else
    templatestrct.specarr{template_ind} = clipstrct.specarr{clip_selectinds};
end


propstrct = clip_propstrct;
set(handles.template_axis,'UserData',propstrct);
set(handles.make_template,'UserData',templatestrct);
plotemplates(handles);

clipstrct.speclabs(clip_selectinds) = templatestrct.speclabs(template_ind);
set(handles.load_clips,'UserData',clipstrct);

plotspecstrct(handles,handles.clip_axis,handles.clip_slider,clipstrct,1000,5,propstrct.srtopt);

deselect_all_Callback(handles.clip_axis);

% --- Executes on button press in remove_template.
function remove_template_Callback(hObject, eventdata, handles)
% hObject    handle to remove_template (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

propstrct = get(handles.template_axis,'UserData');
templatestrct = get(handles.make_template,'UserData');

removeinds = find(propstrct.selectvc);

remlab = templatestrct.speclabs(removeinds);

remtxt = [];
for i = 1:length(removeinds)
    remtxt = [remtxt remlab{i}];
    if i < length(removeinds)
        remtxt = [remtxt ', '];     
    end
end

verif = questdlg(['Verify: about to remove ' remtxt],'Verification','Yes','No','Yes');

if strcmp(verif,'No')
    return;
end

if ~isempty(removeinds)
    templateinds = setdiff(1:length(templatestrct.speclabs),removeinds);
    templatestrct.speclabs = templatestrct.speclabs(templateinds);
    templatestrct.specarr = templatestrct.specarr(templateinds);
    templatestrct.threshvc = templatestrct.threshvc(templateinds);
end

set(handles.make_template,'UserData',templatestrct);

propstrct = get(handles.template_axis,'UserData');
propstrct.selectvc = zeros(1,length(templatestrct.specarr));

set(handles.template_axis,'UserData',propstrct);
plotemplates(handles)

clipstrct = get(handles.load_clips,'UserData');

if ~isempty(clipstrct)
    matchinds = [];
    for remind = 1:length(removeinds)
        matchinds = [matchinds;find(strcmp(clipstrct.speclabs,remlab{remind}))'];
    end

    clipstrct.speclabs(matchinds) = {'x'};
    set(handles.load_clips,'UserData',clipstrct);
    plotspecstrct(handles,handles.clip_axis,handles.clip_slider,clipstrct,1000,5,propstrct.srtopt);

end

% --- Executes on button press in make_class.
function add_to_class_Callback(hObject, eventdata, handles)
% hObject    handle to make_class (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

template_propstrct = get(handles.template_axis,'UserData');
templatestrct = get(handles.make_template,'UserData');
clipstrct = get(handles.load_clips,'UserData');
clip_propstrct = get(handles.clip_axis,'UserData');

clip_selectinds = find(clip_propstrct.selectvc);
template_selectind = find(template_propstrct.selectvc);

if length(template_selectind) ~= 1
    errordlg('Must choose one class');
    return;
end

clipstrct.speclabs(clip_selectinds) = templatestrct.speclabs(template_selectind);

set(handles.load_clips,'UserData',clipstrct);

plotspecstrct(handles,handles.clip_axis,handles.clip_slider,clipstrct,1000,5,clip_propstrct.srtopt);
deselect_all_Callback(handles.clip_axis);



% --- Executes on button press in remove_from_class.
function remove_from_class_Callback(hObject, eventdata, handles)
% hObject    handle to remove_from_class (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

clipstrct = get(handles.load_clips,'UserData');
clip_propstrct = get(handles.clip_axis,'UserData');

clip_selectinds = find(clip_propstrct.selectvc);

clipstrct.speclabs(clip_selectinds) = {'x'};

set(handles.load_clips,'UserData',clipstrct);

plotspecstrct(handles,handles.clip_axis,handles.clip_slider,clipstrct,1000,5,clip_propstrct.srtopt);
deselect_all_Callback(handles.clip_axis);


% --- Executes on button press in sort_len.
function sort_len_Callback(hObject, eventdata, handles)
% hObject    handle to sort_len (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

clipstrct = get(handles.load_clips,'UserData');
plotspecstrct(handles,handles.clip_axis,handles.clip_slider,clipstrct,1000,5,'length');


% --- Executes on button press in sort_lab.
function sort_lab_Callback(hObject, eventdata, handles)
% hObject    handle to sort_lab (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

clipstrct = get(handles.load_clips,'UserData');
plotspecstrct(handles,handles.clip_axis,handles.clip_slider,clipstrct,1000,5,'label');



% --- Executes on button press in load_set.
function load_set_Callback(hObject, eventdata, handles)
% hObject    handle to load_set (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[filename, pathname] = uigetfile('*.mat');
if ~filename
    return;
else
    load([pathname filename])
end

if ~exist(clipstrct.wavdir)
    clipstrct.wavdir = uigetdir(pwd, 'Cannot find original directory');
    
    saveopt = questdlg('Save new directory information?','Save prompt','Yes','No','Yes');
    if strcmp(saveopt,'Yes')
        save([pathname filename],'clipstrct','propstrct','templatestrct')
    end
end

flg = check_headers(clipstrct.wavdir);
clipstrct2 = clipstrct;

load([clipstrct2.wavdir filesep 'clipdirinfo.mat'])
load([clipstrct2.wavdir filesep 'wavdirinfo.mat'])

wavnm = length(dirstrct.wavinds);
clipnm_tot = length(clipstrct.clipsamps);

clear dirstrct

clipstrct = clipstrct2;

clipstrct.clipfl = [pathname filename];

freqs = propstrct.fs*[0:propstrct.f_winlen/2]/propstrct.f_winlen;
if ~isfield(clipstrct,'featarr')
    for clipind = 1:length(clipstrct.specarr)
        %         clipstrct.featarr{clipind} = specfeat(clipstrct.specarr{clipind},freqs(clipstrct.freqinds));
        cliptmp = tdft(clipstrct.wavarr{clipind},propstrct.f_winlen,propstrct.f_winadv);
        cliptmp = cliptmp(clipstrct.freqinds,:);
        clipstrct.featarr{clipind} = specfeat(cliptmp,freqs(clipstrct.freqinds));
    end
end

if ~isfield(clipstrct,'distmat')
    distopt = questdlg('Compute clip distance matrix?','Distance prompt','Yes','No','Yes');
    if strcmp(distopt,'Yes')
        for clipind = 1:length(clipstrct.specarr)
            clipstrct.specarr2{clipind} = tdft(clipstrct.wavarr{clipind},propstrct.f_winlen,propstrct.f_winlen);
            clipstrct.specarr2{clipind} = log(max(clipstrct.specarr2{clipind}(clipstrct.freqinds,:),propstrct.amp_cutoff)) - log(propstrct.amp_cutoff);
        end
        clipstrct.distmat = calcdistmat(clipstrct.specarr2);
        clipstrct.specarr2 = {};
    end
end

set(handles.clip_axis,'UserData',propstrct);
set(handles.load_clips,'UserData',clipstrct);

plotspecstrct(handles,handles.clip_axis,handles.clip_slider,clipstrct,1000,5,'label');

propstrct = get(handles.clip_axis,'UserData');

template_xscope = 1000 / (propstrct.f_winadv * 1000 / propstrct.fs);

set(handles.clip_slider,'Value',1);
clip_slider_Callback(handles.clip_slider, 0, handles);

clipstrct.clipnm_tot = clipnm_tot;
clipstrct.wavnm = wavnm;


titletxt = [pathname filename ' (' num2str(clipstrct.wavnm) ' wav files, ' ...
    num2str(clipstrct.clipnm_tot) ' total clips, ' num2str(length(clipstrct.speclabs)) ' samples)'];

set(handles.title_bar,'String',titletxt);

if isfield(clipstrct,'templatefl') && ~exist('templatestrct')
    load(clipstrct.templatefl)
    templatestrct.speclabs = templatenms;
    templatestrct.threshvc = threshvc;
    templatestrct.templatefl = clipstrct.templatefl;
    templatestrct.wavdir = wavdir;
    templatestrct.freqinds = freqinds;
    
    for tempind = 1:length(templatenms)
        eval(['templatestrct.specarr{tempind} = template' templatenms{tempind} ';'])
    end
    
    if exist('errvc')
        templatestrct.errvc = errvc;
        templatestrct.falsenegvc = falsenegvc;
        templatestrct.temparr = temparr;
        templatestrct.posnonmax = posnonmax;
    end
    
end

set(handles.make_template,'UserData',templatestrct)

propstrct.selectvc = zeros(1,length(templatestrct.speclabs));
propstrct = rmfield(propstrct, 'h');
set(handles.template_axis,'UserData',propstrct);

plotemplates(handles)


% --- Executes on button press in play_clip.
function play_clip_Callback(hObject, eventdata, handles)
% hObject    handle to play_clip (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


slow_fact = str2double(get(handles.slow_factor,'String'));
propstrct = get(handles.clip_axis,'UserData');
selectinds = find(propstrct.selectvc);

if isempty(selectinds)
   errordlg('Must select at least one')
   return;
end

clipstrct = get(handles.load_clips,'UserData');
fs = clipstrct.fs;
selectnm = length(selectinds);

if selectnm == 1
    s = clipstrct.wavarr{selectinds};
    
    if slow_fact ~= 1
         s = pvoc(s, 1/slow_fact, propstrct.f_winlen);
    end
    
else
        
    pad = floor(fs/4);
    wavarrtmp = clipstrct.wavarr(selectinds);
    
    if slow_fact == 1
        slenvc = clipstrct.clipsamps(selectinds);
    else
        slenvc = zeros(1,selectnm);
        
        for selectind = 1:selectnm
            wavarrtmp{selectind} = pvoc(wavarrtmp{selectind}, 1/slow_fact, propstrct.f_winlen);
            slenvc(selectind) = length(wavarrtmp{selectind});           
        end
        
    end
    
    slen = sum(slenvc) + pad*(selectnm-1);
    
    s = zeros(1,slen);
    k = 1;
    
    for selectind = 1:selectnm
        s(k:k+slenvc(selectind)-1) = wavarrtmp{selectind};
        k = k + slenvc(selectind) + pad;
    end
    
end

s = s - mean(s);
slim = max(abs(s));
s = s / slim;

sound(s,fs);


%--------------------------------------------------------------------------
% --- Executes on button press in play_clip.
function extract_clip_Callback(hObject, eventdata, handles)

[filename, pathname] = uiputfile('*.mat');
if ~filename
    return;
else
    fl = [pathname filename];
end

clipstrct = get(handles.load_clips,'UserData');
clip_propstrct = get(handles.clip_axis,'UserData');

clip_selectinds = find(clip_propstrct.selectvc);

spec = clipstrct.specarr{clip_selectinds};
save(fl,'spec');


%--------------------------------------------------------------------------
function change_label_Callback(hObject, eventdata, handles)

propstrct = get(handles.template_axis,'UserData');
templatestrct = get(handles.make_template,'UserData');
clipstrct = get(handles.load_clips,'UserData');

chngind = find(propstrct.selectvc);

if length(chngind) ~= 1
    errordlg('Choose one')
end

classnm = inputdlg('Class label');
if ~isempty(classnm)
    
    if sum(strcmp(templatestrct.speclabs,classnm))>0
        errordlg('Label already in use, please choose another');
        return;
    end
    
    clipinds = find(strcmp(clipstrct.speclabs,templatestrct.speclabs(chngind)));
    templatestrct.speclabs(chngind) = classnm;
    clipstrct.speclabs(clipinds) = classnm;
    
    set(handles.make_template,'UserData',templatestrct);
    set(handles.load_clips,'UserData',clipstrct);
    
    deselect_all_Callback(handles.template_axis);
    plotemplates(handles)
    plotspecstrct(handles,handles.clip_axis,handles.clip_slider,clipstrct,1000,5,propstrct.srtopt);
    
end

%--------------------------------------------------------------------------
function deselect_all_clips_Callback(hObject, eventdata, handles)

deselect_all_Callback(handles.clip_axis);




%--------------------------------------------------------------------------
function import_templates_Callback(hObject, eventdata, handles)
% hObject    handle to load_templates (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[filename, pathname] = uigetfile('*.mat');
if ~filename
    return;
else
    load([pathname filename],'templatestrct')
end

clipstrct = get(handles.load_clips,'UserData');

set(handles.make_template,'UserData',templatestrct)

propstrct = get(handles.clip_axis,'UserData');
propstrct = rmfield(propstrct, 'h');
propstrct.selectvc = zeros(1,length(templatestrct.speclabs));
set(handles.template_axis,'UserData',propstrct);

plotemplates(handles)



% --- Executes on button press in remove_template.
function merge_classes_Callback(hObject, eventdata, handles)
% hObject    handle to remove_template (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

propstrct = get(handles.template_axis,'UserData');
templatestrct = get(handles.make_template,'UserData');

mergeinds = find(propstrct.selectvc);

if length(mergeinds) < 2
    errordlg('Select at least 2');
    return;
end

mergelabs = templatestrct.speclabs(mergeinds);

mergetxt = [];
for i = 1:length(mergeinds)
    mergetxt = [mergetxt mergelabs{i}];
    if i < length(mergeinds)
        mergetxt = [mergetxt ', '];     
    end
end

verif = questdlg(['Verify: about to merge ' mergetxt],'Verification','Yes','No','Yes');

if strcmp(verif,'No')
    return;
else
   warndlg('Any optimization will need to be redone'); 
end

mergelab = inputdlg('Pick new label for merged class (must be one of selected)');
keepind = find(strcmp(mergelabs,mergelab));
if isempty(keepind)
    errordlg('Did not choose one of selected');
    return;
end

set(handles.make_template,'UserData',templatestrct);

propstrct = get(handles.template_axis,'UserData');
propstrct.selectvc = zeros(1,length(templatestrct.specarr));

set(handles.template_axis,'UserData',propstrct);
plotemplates(handles)

clipstrct = get(handles.load_clips,'UserData');

if ~isempty(clipstrct)
    matchinds = [];
    for mergeind = 1:length(mergeinds)
        matchinds = [matchinds;find(strcmp(clipstrct.speclabs,mergelabs{mergeind}))'];
    end

    clipstrct.speclabs(matchinds) = mergelab;
    set(handles.load_clips,'UserData',clipstrct);
    plotspecstrct(handles,handles.clip_axis,handles.clip_slider,clipstrct,1000,5,propstrct.srtopt);

end

templateinds = setdiff(1:length(templatestrct.speclabs),setdiff(mergeinds,mergeinds(keepind)));
templatestrct.speclabs = templatestrct.speclabs(templateinds);
templatestrct.specarr = templatestrct.specarr(templateinds);
templatestrct.threshvc = templatestrct.threshvc(templateinds);

set(handles.make_template,'UserData',templatestrct);

propstrct = get(handles.template_axis,'UserData');
propstrct.selectvc = zeros(1,length(templatestrct.specarr));

set(handles.template_axis,'UserData',propstrct);
plotemplates(handles)




% --- Executes on button press in remove_template.
function remove_from_sample_Callback(hObject, eventdata, handles)
% hObject    handle to remove_template (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

clipstrct = get(handles.load_clips,'UserData');
clip_propstrct = get(handles.clip_axis,'UserData');


verif = questdlg(['Verify: about to remove clips from sample'],'Verification','Yes','No','Yes');

if strcmp(verif,'No')
    return;
end

clip_selectinds = find(clip_propstrct.selectvc);
indsremain = setdiff(1:length(clipstrct.speclabs),clip_selectinds);

clipstrct.speclabs = clipstrct.speclabs(indsremain);
clipstrct.specarr = clipstrct.specarr(indsremain);

if isfield(clipstrct,'specarr2') & ~isempty(clipstrct.specarr2)
    clipstrct.specarr2 = clipstrct.specarr2(indsremain);
end

clipstrct.wavarr = clipstrct.wavarr(indsremain);
if isfield(clipstrct,'distmat')
    clipstrct.distmat = clipstrct.distmat(indsremain,indsremain);
end
% 
% clipstrct.featarr = clipstrct.featarr(indsremain);

clipstrct.clipsamps = clipstrct.clipsamps(indsremain);
clipstrct.clipons = clipstrct.clipons(indsremain);
clipstrct.clipinds = clipstrct.clipinds(indsremain);
clipstrct.wavinds = clipstrct.wavinds(indsremain);
clipstrct.speclens = clipstrct.speclens(indsremain);


clip_propstrct.selectvc = clip_propstrct.selectvc(indsremain);
clip_propstrct.ytops = clip_propstrct.ytops(indsremain);
clip_propstrct.xpos = clip_propstrct.xpos(indsremain);
clip_propstrct.h = clip_propstrct.h(indsremain);
clip_propstrct.speclabs = clip_propstrct.speclabs(indsremain);

set(handles.load_clips,'UserData',clipstrct);
set(handles.clip_axis,'UserData',clip_propstrct);

plotspecstrct(handles,handles.clip_axis,handles.clip_slider,clipstrct,1000,5,clip_propstrct.srtopt);