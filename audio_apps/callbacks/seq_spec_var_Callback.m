function seq_spec_var_Callback(hObject, eventdata, handles)


seqstrct = get(handles.compile_sequences,'UserData');
propstrct = get(get(handles.analysis_panel,'Children'),'UserData');
matchstrct = get(handles.load_match_data,'UserData');

templatestrct = matchstrct.templatestrct;
spec_propstrct = get(handles.template_axis,'UserData');

N = str2num(get(handles.N,'String'));

load([matchstrct.wavdir filesep 'wavdirinfo.mat'])

seqind = find(propstrct.selectvc);
seq = seqstrct.seqs{seqind};
seqvc = zeros(1,length(seq));
seqvc_template = seqvc;

for i = 1:length(seq)
    featms{i} = templatestrct.featms{find(strcmp(seq{i},templatestrct.speclabs))};
    seqvc(i) = find(strcmp(seq{i},seqstrct.labsu));
    seqvc_template(i) = find(strcmp(seq{i},templatestrct.speclabs));
end

[distmat,spmnarr,spstdarr,featms,seqarr,sampinds,spalgnarr] = get_seq_specvar2(seq,seqvc,seqvc_template,featms,seqstrct,matchstrct.wavdir,spec_propstrct,templatestrct,dirstrct,N);
specstrct.wavinds = seqstrct.wavinds(findSeq(seqstrct.vc,seqvc));
specstrct.wavinds = specstrct.wavinds(sampinds);
specstrct.distmat = distmat;
specstrct.seq = seq;
specstrct.spmnarr = spmnarr;
specstrct.spstdarr = spstdarr;
specstrct.cv = zeros(size(seq));
specstrct.std = zeros(size(seq));
specstrct.mn = zeros(size(seq));
specstrct.featms = featms;
specstrct.seqarr = seqarr;
specstrct.sampinds = sampinds;
specstrct.spalgnarr = spalgnarr;

for i = 1:length(seq)
   specstrct.cv(i) = sqrt(mean(vec(specstrct.spstdarr{i}.^2))) /  mean(vec(specstrct.spmnarr{i})); 
   specstrct.cv2(i) = sqrt(sum(vec(specstrct.spstdarr{i}.^2)) / sum(vec(specstrct.spmnarr{i}.^2))); 
   specstrct.std(i) = sqrt(mean(vec(specstrct.spstdarr{i}.^2))); 
   specstrct.mn(i) = mean(vec(specstrct.spmnarr{i})); 
end


set(handles.spec_dist,'UserData',specstrct);
