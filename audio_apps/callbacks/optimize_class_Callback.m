function optimize_class_Callback(hObject, eventdata, handles)
% hObject    handle to make_class (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

template_propstrct = get(handles.template_axis,'UserData');
templatestrct = get(handles.make_template,'UserData');
clipstrct = get(handles.load_clips,'UserData');
clip_propstrct = get(handles.clip_axis,'UserData');

templateind = find(template_propstrct.selectvc);
clipinds = find(strcmp(clipstrct.speclabs,templatestrct.speclabs(templateind)));
clipnm = length(clipinds);
scrs = zeros(clipnm,1);

mininds = 1:clipnm;
mininds2 = mininds;
metaind = 1;
cont = 1;

while cont == 1

    if length(mininds) > 1
        temp1{metaind} = mkTemplateDTW(clipstrct.specarr(clipinds(mininds)));
    else
        temp1(metaind) = clipstrct.specarr(clipinds(mininds));
    end

    for clipind = 1:clipnm
        clip = clipstrct.specarr{clipinds(clipind)};

        scrs(clipind,metaind) = dtwscore(temp1{metaind},clip);
    end

    thresh = prctile(scrs(:,metaind),25);
    mininds2 = find(scrs(:,metaind) < thresh);
    
    if metaind > 1
        mininds = intersect(mininds,find(scrs(:,metaind)<thresh));
    else
        mininds = mininds2;
    end
    %
    %     mininds = find(scrs(:,metaind) < thresh);

    scrs2 = zeros(clipnm,metaind);
    scrs2(:,1:metaind) = scrs;
    scrs = scrs2;
    %
    %     mininds2 = intersect(mininds2,mininds);

    if ~isempty(mininds)
        metaind = metaind + 1;
    else
        cont = 0;
    end

end

templatestrct.temparr{templateind} = temp1;
set(handles.make_template,'UserData',templatestrct);
deselect_all_Callback(handles.template_axis);