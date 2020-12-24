function plot_class_specdist_Callback(hObject, eventdata, handles)
% hObject    handle to plot_tempmatches (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

template_propstrct = get(handles.template_axis,'UserData');
templatestrct = get(handles.make_template,'UserData');
clipstrct = get(handles.load_clips,'UserData');
clip_propstrct = get(handles.clip_axis,'UserData');

template_selectinds = find(template_propstrct.selectvc);

h = subplot(1,1,1,'Parent',handles.analysis_panel);
cla(h,'reset');

if ~isfield(clipstrct,'distmat')
    distopt = questdlg('Must calculate distance matrix first, do it now?','Distance prompt','Yes','No','Yes')
    
    if strcmp(distopt,'Yes')
        for clipind = 1:length(clipstrct.specarr)
            clipstrct.specarr2{clipind} = tdft(clipstrct.wavarr{clipind},clip_propstrct.f_winlen,clip_propstrct.f_winlen);
            clipstrct.specarr2{clipind} = log(max(clipstrct.specarr2{clipind}(clipstrct.freqinds,:),clip_propstrct.amp_cutoff)) - log(clip_propstrct.amp_cutoff);
        end
        clipstrct.distmat = calcdistmat(clipstrct.specarr2);
        
        set(handles.load_clips,'UserData',clipstrct)
        
    else  
        return;
    end
end


labstmp = clipstrct.speclabs;
indsval = [];

for ind = 1:length(template_selectinds)
    indstmp = find(strcmp(labstmp,templatestrct.speclabs(template_selectinds(ind))));
    indsval = union(indsval,indstmp);
end

indstmp = find(strcmp(labstmp,'x'));
indsval = union(indsval,indstmp);

labstmp = labstmp(indsval);

labarr = unique(labstmp);
labvc = zeros(length(indsval),1);
for labind = 1:length(labarr)
    labvc(find(strcmp(labstmp,labarr{labind}))) = labind;
end

indsmat = repmat(labvc,1,length(labstmp));
indsmat2 = repmat([1:length(labstmp)]',1,length(labstmp));

dimnm = length(labarr);
distmat = clipstrct.distmat(indsval,indsval);

h = zeros(dimnm);
xlims = [0 1.1*max(distmat(find(indsmat2 > indsmat2')))];

xvc = 0:.5:xlims(end);

for dimind1 = 1:dimnm
    for dimind2 = 1:dimnm
     
        if dimind1==dimind2
            clr = 'b';
        else
            clr = 'r';
        end
        
        h(dimind2,dimind1)=subplot(dimnm,dimnm,dimind1+(dimind2-1)*dimnm,'Parent',handles.analysis_panel);
        indstmp = find(indsmat2 ~= indsmat2' & indsmat==dimind1 & indsmat'==dimind2);
        [N,xout] = hist(distmat(indstmp),xvc);
        bar(h(dimind2,dimind1),xout,N/max(N),'linestyle','none','facecolor',clr);
        
        xlim(xlims);
        
        if dimind2 < dimnm
            set(h(dimind2,dimind1),'xtick',[])
        end
        
        if dimind1 > 1
            set(h(dimind2,dimind1),'ytick',[])
        end
        
        if dimind2 == dimnm
            xlabel(labarr{dimind1});
        end
        
        if dimind1 == 1
            ylabel(labarr{dimind2})
        end
    
        
        ylim([0 1.1])
        set(h(dimind2,dimind1),'ytick',[],'box','off');
        
        
    end
    
end

if dimnm == 1
   ylabel('count')
%    set(gca,'ytick',[0:.1:1],'yticklabel',[0:.1:1])
    
end
