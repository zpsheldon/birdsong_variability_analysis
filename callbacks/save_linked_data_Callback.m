function save_linked_data_Callback(hObject, eventdata, handles)

matchstrct = get(handles.load_match_data,'UserData');
templatestrct = matchstrct.templatestrct;

matchstrct = rmfield(matchstrct,'templatestrct');
matchstrct = rmfield(matchstrct,'trainingfl');

[matchfl,matchpth] = uiputfile('*.mat','Linked match file name');
[templatefl,templatepth] = uiputfile('*.mat','Linked template file name');

matchstrct.templatefl = [templatepth templatefl];

save([matchpth matchfl],'matchstrct');
templatenms = templatestrct.speclabs(1:end-1);

propstrct = get(handles.template_axis,'UserData');
winlen = propstrct.f_winlen;
winadv = propstrct.f_winadv;
filt_f1 = propstrct.freqmin;
filt_f2 = propstrct.freqmax;
filt_order = propstrct.filt_order;
fs = propstrct.fs;
freqinds = templatestrct.freqinds;
threshvc = templatestrct.threshvc;
wavdir = templatestrct.wavdir;

save([templatepth templatefl],...
    'winlen','winadv','filt_f1','filt_f2','filt_order','fs','freqinds','threshvc','templatenms','wavdir')

for tempind = 1:length(templatenms)
    eval(['template' templatenms{tempind} '=templatestrct.specarr{' num2str(tempind) '};'])
    save([templatepth templatefl],['template' templatenms{tempind}],'-append');
end
