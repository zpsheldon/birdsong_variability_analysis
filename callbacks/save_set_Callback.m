function save_set_Callback(hObject, eventdata, handles)
% hObject    handle to save_set (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

propstrct = get(handles.clip_axis,'UserData');
clipstrct = get(handles.load_clips,'UserData');
templatestrct = get(handles.make_template,'UserData');

if isfield(clipstrct,'clipfl')
    try
        save(clipstrct.clipfl,'propstrct','clipstrct','templatestrct');        
    catch
        warndlg('Did not save, permission denied by Matlab to overwrite (try changing workspace directory')
    end
else
    saveas_set_Callback(0, 0, handles)
end