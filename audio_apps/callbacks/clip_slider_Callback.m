function clip_slider_Callback(hObject, eventdata, handles)
% hObject    handle to clip_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

ypos = get(hObject, 'Value');
propstrct = get(handles.clip_axis,'UserData');
ymin = propstrct.ymin + propstrct.ymax;
ymax = propstrct.ymax;
valmax = ymin + ypos*(ymax - ymin);
valmin = valmax - propstrct.yscope;

axes(handles.clip_axis)

inds = find(propstrct.ytops >= valmin & propstrct.ytops <= valmax & ~propstrct.labeled);
for clipind = 1:length(inds)
    text_x = propstrct.xpos(inds(clipind));
    text_y = propstrct.ytops(inds(clipind));
    text(text_x,text_y,propstrct.speclabs(inds(clipind)),'HorizontalAlignment','center','VerticalAlignment','bottom','Clipping','on');
end
propstrct.labeled(inds) = 1;

set(handles.clip_axis,'ylim',[valmin valmax],'UserData',propstrct)
refresh(handles.mainfig);