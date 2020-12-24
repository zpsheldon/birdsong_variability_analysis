function deselect_all_Callback(select_axis)
% hObject    handle to deselect_all (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

propstrct = get(select_axis,'UserData');
selectinds = find(propstrct.selectvc);

for i = 1:length(selectinds)
%     set(propstrct.h(selectinds(i)),'Selected','off');
    clipick(propstrct.h(selectinds(i)));
end

propstrct.selectvc(selectinds) = 0;
set(select_axis,'UserData',propstrct);
