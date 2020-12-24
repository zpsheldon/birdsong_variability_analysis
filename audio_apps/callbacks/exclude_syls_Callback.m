function exclude_syls_Callback(hObject, eventdata, handles)

propstrct = get(handles.template_axis,'UserData');
template_selectinds = find(propstrct.selectvc);

matchstrct = get(handles.load_match_data,'UserData');
templatestrct = matchstrct.templatestrct;


sylset = get(handles.exclude_syls,'UserData');

sylset(template_selectinds) = 1 - sylset(template_selectinds);

plotspecstrct(handles,handles.template_axis,[],templatestrct,1000,0,'label');
propstrct = get(handles.template_axis,'UserData');


for ind = 1:length(sylset)
    if ~sylset(ind)
        CData = get(propstrct.h(ind),'CData');
        XData = get(propstrct.h(ind),'XData');
        YData = get(propstrct.h(ind),'YData');
        CData = 2+CData*.15;
        XData = median(XData) + (XData-median(XData))*.8;
        YData = max(YData) + (YData-max(YData))*.9;
        set(propstrct.h(ind),'CData',CData,'XData',XData,'YData',YData);
    end
end

set(handles.exclude_syls,'UserData',sylset);

deselect_all_Callback(handles.template_axis);

calc_seqvc(handles);

transition_matrix_Callback(0,0,handles);