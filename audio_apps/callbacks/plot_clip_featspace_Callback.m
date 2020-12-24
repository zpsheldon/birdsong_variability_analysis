function plot_clip_featspace_Callback(hObject, eventdata, handles)
% hObject    handle to plot_tempmatches (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

template_propstrct = get(handles.template_axis,'UserData');
templatestrct = get(handles.make_template,'UserData');
clipstrct = get(handles.load_clips,'UserData');
clip_propstrct = get(handles.clip_axis,'UserData');

template_selectinds = find(template_propstrct.selectvc);
clip_selectinds = find(clip_propstrct.selectvc);

h = subplot(1,1,1,'Parent',handles.analysis_panel);
cla(h,'reset');

labarr = ['x' templatestrct.speclabs(template_selectinds)];
labnm = length(labarr);

clipnm = length(clipstrct.speclabs);
featnm = size(clipstrct.featarr{1},1)+1;
featmat = zeros(clipnm,featnm);

featlabs = {'length (msec)','avg fund freq (Hz)','avg ceps Wiener ent','avg log-amp'};

dt = clip_propstrct.f_winlen * 1000 / clip_propstrct.fs;

for clipind = 1:clipnm
    featmat(clipind,:) = [clipstrct.speclens(clipind)*dt mean(clipstrct.featarr{clipind}')];
end

labstmp = clipstrct.speclabs;
indsnon = 1:length(labstmp);
labarr(end+1) = {'clips'};

for ind = 1:length(template_selectinds)
    indstmp = find(~strcmp(labstmp,templatestrct.speclabs(template_selectinds(ind))));
    indsnon = intersect(indsnon,indstmp);
end

labstmp(indsnon) = {'x'};

markarr = {'o','s','d'};
clrarr = {'r','b','g','c','m','y'};

marknm = length(markarr);
clrnm = length(clrarr);

symnm = marknm * clrnm;

clrarr = ['k' repmat(clrarr,1,marknm)];
markarr = ['.' repmat(markarr,1,clrnm)];

dimnm = featnm;

h = zeros(dimnm);
h_pts = zeros(size(featmat));

for dimind1 = 1:dimnm
    for dimind2 = 1:dimnm
        
        h(dimind2,dimind1)=subplot(dimnm,dimnm,dimind1+(dimind2-1)*dimnm,'Parent',handles.analysis_panel);
        
        bufftmp = .05*abs(max(featmat(:,dimind1)) - min(featmat(:,dimind1)));
        xlims = [min(featmat(:,dimind1))-bufftmp max(featmat(:,dimind1))+bufftmp];
        
        hold on
        
        if dimind1 ~= dimind2           
          
            bufftmp = .05*abs(max(featmat(:,dimind2)) - min(featmat(:,dimind2)));
            ylims = [min(featmat(:,dimind2))-bufftmp max(featmat(:,dimind2))+bufftmp];
            ylim(ylims);
           
            valmax = min(ylims(end),xlims(end));
              
            for labind = 1:labnm            
                indstmp = find(strcmp(labstmp,labarr(labind)));  
                x1 = featmat(indstmp,dimind1);
                x2 = featmat(indstmp,dimind2);
                htmp = plot(h(dimind2,dimind1),x1,x2,...
                    'Marker',markarr{labind},'MarkerFaceColor','none','MarkerEdgeColor',clrarr{labind},'LineStyle','none','MarkerSize',4);
            end
  
                        
            plot(featmat(clip_selectinds,dimind1),featmat(clip_selectinds,dimind2),'ko','MarkerSize',6,'MarkerFaceColor','g');
            

        else
          
            for labind = 1:labnm
                indstmp = find(strcmp(labstmp,labarr(labind)));
                x = featmat(indstmp,dimind1);
                [N,xout] = hist(x);
                plot(xout,N/max(N),[markarr{labind} '-' clrarr{labind}],'linewidth',1.5,'MarkerSize',4)
                hold on
            end
            
            plot([featmat(clip_selectinds,dimind1) featmat(clip_selectinds,dimind1)],[0 1.1],'g-','LineWidth',2);
            
            
            if dimind1 == 1 && dimind2 == 1
                htmp = legend(labarr,'Location','Best','fontsize',6);
            end
            
            ylim([0 1.1])
                        
        end
        
        xlim(xlims);
        
        if dimind2 < dimnm
            set(h(dimind2,dimind1),'xtick',[])
        end
        
        if dimind1 > 1
            set(h(dimind2,dimind1),'ytick',[])
        end
        
        if dimind2 == dimnm
            xlabel(featlabs{dimind1});
        end
        
        if dimind1 == 1
            ylabel(featlabs{dimind2})
        end
        
        set(h(dimind2,dimind1),'FontSize',8);
        
        hold off
        
    end
    
end

if dimnm > 1
    yticks = get(h(end,1),'xtick');
    fct = diff(get(h(1,1),'ylim')) / diff(get(h(1,end),'ylim'));
    set(h(1,1),'ytick',yticks*fct,'yticklabel',yticks)
end


if dimnm == 1
   ylabel('count')
%    set(gca,'ytick',[0:.1:1],'yticklabel',[0:.1:1])
    
end

try
    [class,err] = classify(featmat,featmat,labstmp','quadratic');
    htmp = get(h(1,2),'Title');
    set(htmp,'String',['quadratic discriminant error = ' num2str(round(1000*err)/10) '%'],'fontweight','bold','Units','normalized','Position',[1 1.3 0],'HorizontalAlignment','center');
catch
end