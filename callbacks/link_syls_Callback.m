function link_syls_Callback(hObject, eventdata, handles)

propstrct = get(handles.template_axis,'UserData');

matchstrct = get(handles.load_match_data,'UserData');
templatestrct = matchstrct.templatestrct;

linkseq = inputdlg('Link sequence (separate syls by commas)');
linklab = inputdlg('Link label');

cliplabs = matchstrct.cliplabs;

sylset = get(handles.exclude_syls,'UserData');
indsval = 1:length(matchstrct.cliplabs);
if ~isempty(sylset)
    exinds = find(sylset==0);
    for exind = 1:length(exinds)
        indsval = intersect(indsval,find(~strcmp(matchstrct.cliplabs,matchstrct.templatestrct.speclabs(exinds(exind)))));
    end
end

[vc,labtmp,seqinds2,clipons2,cliplens2] = labs2vc(matchstrct.cliplabs(indsval),matchstrct.wavinds(indsval),...
    matchstrct.cliptms(indsval),matchstrct.cliplens(indsval),str2double(get(handles.max_gap,'String')));

linkseq = linkseq{1};
commainds = [0 findstr(linkseq,',') length(linkseq)+1];
linkvc = zeros(1,length(commainds)-1);
linkvc2 = linkvc;
for linkind = 1:length(linkvc)
    linkvc(linkind) = find(strcmp(labtmp,linkseq(commainds(linkind)+1:commainds(linkind+1)-1)));
    linkvc2(linkind) = find(strcmp(templatestrct.speclabs,linkseq(commainds(linkind)+1:commainds(linkind+1)-1)));
end


linkinds = findSeq(vc,linkvc);
ons_new = clipons2(linkinds);
offs_new = clipons2(linkinds+length(linkvc)-1) + cliplens2(linkinds+length(linkvc)-1);
wavinds_new = seqinds2(linkinds);

matchstrct_new = matchstrct;
clipoffs = matchstrct_new.cliptms + matchstrct_new.cliplens;
for matchind = 1:length(linkinds)
    
    delinds = find(matchstrct_new.wavinds == wavinds_new(matchind) & matchstrct_new.cliptms >= ons_new(matchind) & clipoffs <= offs_new(matchind));
    delindsval = [];
    for delind = 1:length(delinds)
        if ~any(strcmp(matchstrct_new.cliplabs(delinds(delind)),linkseq))
            delindsval = intersect(delindsval,delind);
        end
    end
 
    matchstrct_new.cliplabs(delinds(1)) = linklab;
    matchstrct_new.cliplens(delinds(1)) = offs_new(matchind) - ons_new(matchind);
    
    indsval = setdiff(1:length(matchstrct_new.wavinds),delinds(2:end));
    
    matchstrct_new.cliplabs = matchstrct_new.cliplabs(indsval);
    matchstrct_new.wavinds = matchstrct_new.wavinds(indsval);
    matchstrct_new.cliplens = matchstrct_new.cliplens(indsval);
    matchstrct_new.cliptms = matchstrct_new.cliptms(indsval);
    matchstrct_new.clipinds = matchstrct_new.clipinds(indsval);
    matchstrct_new.clipnm = matchstrct_new.clipnm - length(delinds) + 1;
    
    clipoffs = matchstrct_new.cliptms + matchstrct_new.cliplens;

end

tempnew = [];

for linkind = 1:length(linkvc2)
    tempnew = [tempnew templatestrct.specarr{linkvc2(linkind)}];
end

%indsnew = setdiff(1:length(templatestrct.speclabs),linkvc2);

indsnew = 1:length(templatestrct.speclabs);
templatestrct.speclabs = [templatestrct.speclabs(indsnew(1:end-1)) linklab 'x'];
templatestrct.specarr = [templatestrct.specarr(indsnew(1:end-1)) tempnew templatestrct.specarr(end)];

plotspecstrct(handles,handles.template_axis,[],templatestrct,1000,0,'label');

matchstrct_new.templatestrct = templatestrct;
set(handles.load_match_data,'UserData',matchstrct_new);

deselect_all_Callback(handles.template_axis);
calc_seqvc(handles);

transition_matrix_Callback(0,0,handles);