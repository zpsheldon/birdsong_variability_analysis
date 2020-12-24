function plotemplates(handles)

propstrct = get(handles.template_axis,'UserData');
templatestrct = get(handles.make_template,'UserData');

propstrct.selectvc = zeros(1,length(templatestrct.specarr));
set(handles.template_axis,'UserData',propstrct);

plotspecstrct(handles,handles.template_axis,[],templatestrct,1000,0,'label');