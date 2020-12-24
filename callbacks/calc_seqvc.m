function calc_seqvc(handles)

matchstrct = get(handles.load_match_data,'UserData');
sylset = get(handles.exclude_syls,'UserData');
indsval = 1:length(matchstrct.cliplabs);
if ~isempty(sylset)
    exinds = find(sylset==0);
    for exind = 1:length(exinds)
        indsval = intersect(indsval,find(~strcmp(matchstrct.cliplabs,matchstrct.templatestrct.speclabs(exinds(exind)))));
    end
end

[vc,labtmp,seqinds2,clipons2,cliplens2] = labs2vc(matchstrct.cliplabs(indsval),matchstrct.wavinds(indsval),matchstrct.cliptms(indsval),matchstrct.cliplens(indsval),...
    str2double(get(handles.max_gap,'String')));

seqstrct.matchinds = indsval;
seqstrct.wavinds = seqinds2;
seqstrct.cliptms = clipons2;
seqstrct.cliplens = cliplens2;
seqstrct.vc = vc;
seqstrct.labsu = labtmp;

set(handles.compile_sequences,'UserData',seqstrct);