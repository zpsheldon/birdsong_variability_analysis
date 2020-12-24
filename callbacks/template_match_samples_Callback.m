function template_match_samples_Callback(hObject, eventdata, handles)
% hObject    handle to match_all_samples (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

clipstrct = get(handles.load_clips,'UserData');
templatestrct = get(handles.make_template,'UserData');
clip_propstrct = get(handles.clip_axis,'UserData');
template_propstrct = get(handles.template_axis,'UserData');

templateinds = find(template_propstrct.selectvc);
templatenm = length(templateinds);

clipinds = find(strcmp(clipstrct.speclabs,'x'));
clipnm = length(clipinds);

matchmat = zeros(clipnm,templatenm);

if isempty(templateinds)
    errordlg('Choose at least 1');
    return;
end

clipstrct.ids = {};
clipstrct.scr = zeros(1,clipnm);
matchtot = clipnm * length(templateinds);
matchind = 1;

h = waitbar(1/matchtot,'Matching samples');

for clipind = 1:clipnm
    clip = clipstrct.specarr{clipinds(clipind)};
    clip = clip / max(clip(:));
    
    for templateindind = 1:templatenm
        
        templateind = templateinds(templateindind);               
        matchmat(clipind,templateindind) = dtwscore(templatestrct.specarr{templateind}/max(templatestrct.specarr{templateind}(:)),clip);         
        
        matchind = matchind + 1;
        h = waitbar(matchind/matchtot,h);

    end

end

close(h)

if templatenm > 1
    [maxscrs,maxinds] = min(matchmat');
else
    maxscrs = matchmat;
    maxinds = ones(clipnm,1);
end

valvc = maxscrs <= .04;
clipstrct.speclabs(clipinds(valvc)) = templatestrct.speclabs(templateinds(maxinds(valvc)));

set(handles.load_clips,'UserData',clipstrct);
plotspecstrct(handles,handles.clip_axis,handles.clip_slider,clipstrct,1000,5,'label');