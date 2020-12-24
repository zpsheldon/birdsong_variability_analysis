function optimize_threshold_Callback(hObject, eventdata, handles)
% hObject    handle to optimize_threshold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

template_propstrct = get(handles.template_axis,'UserData');
templatestrct = get(handles.make_template,'UserData');
clipstrct = get(handles.load_clips,'UserData');
clip_propstrct = get(handles.clip_axis,'UserData');

template_selectinds = find(template_propstrct.selectvc);
[dmy,maxinds] = max(clipstrct.matchmat,[],2);

xvc = 0:.1:3;

distrct = get(handles.sample_clip_matches,'UserData');
if isempty(distrct)
    distvc = ones(1,length(xvc));
else
    distmat = distrct.distmat2;
    distvc = distmat(:,template_selectinds)';
end

for i = 1:length(template_selectinds)

    maxinds2 = find(maxinds==template_selectinds(i));

    matchinds = find(strcmp(clipstrct.speclabs(maxinds2),templatestrct.speclabs(template_selectinds(i))));
    nonmatchinds = setdiff(1:length(maxinds2),matchinds);

    N_match = histc((clipstrct.matchmat(maxinds2(matchinds),template_selectinds(i))),xvc)';
    N_nonmatch = histc((clipstrct.matchmat(maxinds2(nonmatchinds),template_selectinds(i))),xvc)';
    
    N_match = distvc(i,:) .* N_match(:)' / sum(N_match(:));
    N_nonmatch = distvc(i,:) .* N_nonmatch(:)' / sum(N_nonmatch(:));


    cumprob = [fliplr(cumsum([0 fliplr(N_match)]));cumsum([0 N_nonmatch])];
    
    correctdec = sum(cumprob);
    optind = find(correctdec==max(correctdec));

    templatestrct.threshvc(template_selectinds(i)) = min(xvc(optind));

end

set(handles.make_template,'UserData',templatestrct);

analyze_templates_Callback(0,0, handles);
plotemplates(handles)