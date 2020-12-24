function plot_clip_features_Callback(hObject, eventdata, handles)
% hObject    handle to plot_clip_matches (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

h = subplot(1,1,1,'Parent',handles.analysis_panel);

reset(handles.analysis_panel);

template_propstrct = get(handles.template_axis,'UserData');
templatestrct = get(handles.make_template,'UserData');
clipstrct = get(handles.load_clips,'UserData');
clip_propstrct = get(handles.clip_axis,'UserData');

clip_selectinds = find(clip_propstrct.selectvc);
clipnm = length(clip_selectinds);

if clipnm == 0
    errordlg('Must select clip(s)');
    return;
end

featnm = size(clipstrct.featarr{1},1);

maxspeclen = max(clipstrct.speclens(clip_selectinds));

featmat = zeros(clipnm,featnm);
featmat2 = zeros(clipnm,featnm,maxspeclen);

dt = clip_propstrct.f_winadv * 1000 / clip_propstrct.fs;


markarr = {'o','s','d'};
clrarr = {'b','r','g','c','m','y'};

marknm = length(markarr);
clrnm = length(clrarr);

symnm = marknm * clrnm;

clrarr = [repmat(clrarr,1,marknm)];
markarr = [repmat(markarr,1,clrnm)];

h = zeros(featnm,1);

for clipind = 1:clipnm
   featmp = clipstrct.featarr{clip_selectinds(clipind)};
   featmat2(clipind,:,1:size(featmp,2)) = featmp;
   featmat(clipind,:) = median(featmp')';
   
   for featind = 1:featnm    
       h(featind) = subplot(4,3,(featind*3+1):(featind*3+3),'Parent',handles.analysis_panel);
       t = dt*[0:clipstrct.speclens(clip_selectinds(clipind))-1];
       plot(t,featmp(featind,:),[clrarr{clipind} '-'])
       hold on
   end
end

for featind = 1:featnm
    set(h(featind),'xlim',[0 (maxspeclen-1)*dt],'box','off')
end

legend(clipstrct.speclabs(clip_selectinds),'Location','South');

set(get(h(1),'YLabel'),'String','fund. freq. (Hz)')
set(get(h(2),'YLabel'),'String','ceps Wiener ent')
set(get(h(3),'YLabel'),'String','log-amplitude')

set(h(1),'xticklabel',[])
set(h(2),'xticklabel',[])
set(get(h(3),'XLabel'),'String','time (msec)')

htmp = subplot(4,3,1,'Parent',handles.analysis_panel);
bar(htmp,featmat(:,1))
set(htmp,'box','off')
ylabel('median fund freq (Hz)')
xlim([.5 clipnm+.5])
set(htmp,'xtick',1:clipnm,'xticklabel',clipstrct.speclabs(clip_selectinds))

htmp = subplot(4,3,2,'Parent',handles.analysis_panel);
bar(htmp,featmat(:,2))
set(htmp,'box','off')
ylabel('ceps Wiener ent')
xlim([.5 clipnm+.5])
set(htmp,'xtick',1:clipnm,'xticklabel',clipstrct.speclabs(clip_selectinds))

htmp = subplot(4,3,3,'Parent',handles.analysis_panel);
bar(htmp,featmat(:,3))
set(htmp,'box','off')
ylabel('median log-amp')
xlim([.5 clipnm+.5])
set(htmp,'xtick',1:clipnm,'xticklabel',clipstrct.speclabs(clip_selectinds))
