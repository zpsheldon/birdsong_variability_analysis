function saveas_set_Callback(hObject, eventdata, handles)
% hObject    handle to save_set (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

propstrct = get(handles.clip_axis,'UserData');
clipstrct = get(handles.load_clips,'UserData');
templatestrct = get(handles.make_template,'UserData');

[filename, pathname] = uiputfile('*.mat');
if filename
    try
        clipstrct.clipfl = [pathname filename];
        save([pathname,filename],'propstrct','clipstrct','templatestrct');
        set(handles.load_clips,'UserData',clipstrct);
        
        if isfield(clipstrct,'wavnm')
            titletxt = [pathname filename ' (' num2str(clipstrct.wavnm) ' wav files, ' ...
                num2str(clipstrct.clipnm_tot) ' total clips, ' num2str(length(clipstrct.speclabs)) ' samples)'];
        else
            titletxt = [pathname filename ' (' num2str(length(clipstrct.speclabs)) ' samples)'];
        end
        
        set(handles.title_bar,'String',titletxt);
        
    catch
        warndlg('Did not save, permission denied by Matlab to overwrite (try changing workspace directory)')
    end
end