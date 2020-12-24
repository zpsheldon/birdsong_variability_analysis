function match_all_samples_Callback(hObject, eventdata, handles)
% hObject    handle to match_all_samples (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

clipstrct = get(handles.load_clips,'UserData');
templatestrct = get(handles.make_template,'UserData');
clip_propstrct = get(handles.clip_axis,'UserData');

templatenm = length(templatestrct.speclabs);
clipnm = length(clipstrct.speclabs);
clipstrct.matchmat = zeros(clipnm,templatenm);
clipstrct.indmat = clipstrct.matchmat;

clipstrct.ids = {};
clipstrct.scr = zeros(1,clipnm);
matchtot = clipnm * templatenm;
matchind = 1;

h = waitbar(1/matchtot,'Matching samples');

for clipind = 1:clipnm
    clip = clipstrct.specarr{clipind};
    
    for templateind = 1:templatenm
        
        matchind = matchind + 1;
        h = waitbar(matchind/matchtot,h);
        
        tempnm = length(templatestrct.temparr{templateind});
        scrtmp = zeros(tempnm,1);
        
        for tempind = 1:tempnm
            scrtmp(tempind) = dtwscore(templatestrct.temparr{templateind}{tempind},clip);
        end

        [clipstrct.matchmat(clipind,templateind),clipstrct.indmat(clipind,templateind)] = max(scrtmp);

%          clipstrct.matchmat(clipind,templateind) =  dtwscore(templatestrct.specarr{templateind},clip);
%         
    end
    
    [scr,idind] = max(clipstrct.matchmat(clipind,:));
    
    if clipstrct.matchmat(clipind,idind) >= templatestrct.threshvc(idind)
        clipstrct.speclabs(clipind) = templatestrct.speclabs(idind);
        clipstrct.scr(clipind) = scr;
    else
        clipstrct.speclabs{clipind} = 'x';
    end
    
end

close(h)

[dmy,srtind] = sort(clipstrct.speclabs);

set(handles.load_clips,'UserData',clipstrct);
plotspecstrct(handles,handles.clip_axis,handles.clip_slider,clipstrct,1000,5,'label');