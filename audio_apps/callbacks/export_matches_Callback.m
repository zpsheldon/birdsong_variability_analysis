function export_matches_Callback(hObject, eventdata, handles)
% hObject    handle to export_matches (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


propstrct = get(handles.template_axis,'UserData');
templatestrct = get(handles.make_template,'UserData');

load([templatestrct.wavdir filesep 'clipdirinfo.mat']);
load([templatestrct.wavdir filesep 'wavdirinfo.mat']);

cliplens = 1000 * double(clipstrct.clipsamps) / propstrct.fs;
cliptms = 1000 * double(clipstrct.clipons) / propstrct.fs;

[scrs,maxinds] = max(templatestrct.matchmat_all');
threshvc = templatestrct.threshvc(maxinds);
indsval = find(scrs >= threshvc);
indsnon = find(scrs < threshvc);

cliplabs(indsval) = templatestrct.speclabs(maxinds(indsval));
cliplabs(indsnon) = {'x'};

[flnm,flpth] = uiputfile('*.mat','Choose where to save match data');
filename = [flpth flnm];
if ~filename
    return;
end

if ~isfield(templatestrct,'templatefl')
    save_templates_Callback(0, 0, handles);
    templatestrct = get(handles.make_template,'UserData');
end

matchstrct.wavdir = templatestrct.wavdir;
matchstrct.wavinds = clipstrct.wavinds;
matchstrct.clipinds = clipstrct.clipinds;
matchstrct.cliptms = cliptms;
matchstrct.clipscrs = scrs;
matchstrct.cliplabs = cliplabs;
matchstrct.cliplens = cliplens;
matchstrct.templatefl = templatestrct.templatefl;

matchmat = templatestrct.matchmat_all;
indmat = templatestrct.indmat_all;

save(filename,'matchstrct','matchmat','indmat')