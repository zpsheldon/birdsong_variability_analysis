function plot_clip_matches_Callback(hObject, eventdata, handles)
% hObject    handle to plot_clip_matches (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

h = subplot(1,1,1,'Parent',handles.analysis_panel);

reset(handles.analysis_panel);

template_propstrct = get(handles.template_axis,'UserData');
templatestrct = get(handles.make_template,'UserData');
clipstrct = get(handles.load_clips,'UserData');
clip_propstrct = get(handles.clip_axis,'UserData');

clip_selectinds = find(clip_propstrct.selectvc);
clipnm = length(clip_selectinds);

tempnm = length(templatestrct.specarr);

if clipnm == 0
    errordlg('Must select clip(s)');
    return;
end

matchmat = zeros(clipnm,tempnm);

for tempind = 1:tempnm
    for clipind = 1:clipnm
        matchmat(tempind,clipind) = ...
            1 / dtwscore(templatestrct.specarr{tempind}/max(templatestrct.specarr{tempind}(:)),clipstrct.specarr{clip_selectinds(clipind)}/max(clipstrct.specarr{clip_selectinds(clipind)}(:)));
    end
end

bar(h,matchmat)

ymax = max(max(matchmat));

set(h,'xtick',1:length(templatestrct.speclabs),'xticklabel',templatestrct.speclabs,...
    'ylim',[0 1.1*ymax],'xlim',[.3 length(templatestrct.speclabs)+.7],'box','off');

if length(clip_selectinds) > 1
    legend(clip_propstrct.speclabs(clip_selectinds),'Location','North','box','off');
end

