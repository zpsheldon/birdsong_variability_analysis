function set_threshold_Callback(hObject, eventdata, handles)
% hObject    handle to set_threshold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

templatestrct = get(handles.make_template,'UserData');
propstrct = get(handles.template_axis,'UserData');

tempind = find(propstrct.selectvc);

thresh = inputdlg('put in threshold');
thresh = str2num(thresh{1});
templatestrct.threshvc(tempind) = thresh;

set(handles.make_template,'UserData',templatestrct);
plotemplates(handles)
