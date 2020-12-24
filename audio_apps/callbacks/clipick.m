function varargout = clipick(hObject,eventdata,handles)

state = get(hObject,'Selected');
propstrct = get(get(hObject,'Parent'),'UserData');
id = get(hObject,'UserData');

cflag = 0;
if isfield(get(hObject),'CData')
    hax = get(hObject,'Parent');
    clims = get(hax,'clim');
    cflag = 1;
end

switch state
    case 'off'
        set(hObject,'Selected','on');
        propstrct.selectvc(id) = 1;
        
        if cflag
            set(hObject,'CData',get(hObject,'CData')+clims(2)/3)
        end

    case 'on'
        set(hObject,'Selected','off');
        propstrct.selectvc(id) = 0;
        
        if cflag
            set(hObject,'CData',get(hObject,'CData')-clims(2)/3)
        end

end

set(get(hObject,'Parent'),'UserData',propstrct);

