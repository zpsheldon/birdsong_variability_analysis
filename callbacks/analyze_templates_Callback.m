function analyze_templates_Callback(hObject, eventdata, handles)
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
xvc = 0:5:100;

labarr = ['x' templatestrct.speclabs(template_selectinds)];
labnm = length(labarr);
threshvc = templatestrct.threshvc(template_selectinds);
matchtmp = clipstrct.matchmat(:,template_selectinds);

[dmy,maxinds] = max(clipstrct.matchmat');

labstmp = clipstrct.speclabs;
indsnon = 1:length(labstmp);

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

dimnm = length(template_selectinds);

h = zeros(dimnm);
h_pts = zeros(size(matchtmp));

for dimind1 = 1:dimnm
    for dimind2 = 1:dimnm
        
        h(dimind2,dimind1)=subplot(dimnm,dimnm,dimind1+(dimind2-1)*dimnm,'Parent',handles.analysis_panel);
        
        hold on
        
        xlims = [0 max(matchtmp(:,dimind1))+.05];
        if dimind1 ~= dimind2           
          
            ylims = [min(matchtmp(:,dimind2))-.05 max(matchtmp(:,dimind2))+.05];
            ylim(ylims);
            
            threshmin = min(threshvc([dimind1 dimind2]));
            valmax = min(ylims(end),xlims(end));
              
            for labind = 1:labnm            
                indstmp = find(strcmp(labstmp,labarr(labind)));  
                htmp = plot(h(dimind2,dimind1),matchtmp(indstmp,dimind1),matchtmp(indstmp,dimind2),...
                    'Marker',markarr{labind},'MarkerFaceColor',clrarr{labind},'MarkerEdgeColor',clrarr{labind},'LineStyle','none');
            end
            
            if dimind1 == 2 & dimind2 == 1
                legend(labarr)
            end
            
            plot(h(dimind2,dimind1),[threshmin valmax],[threshmin,valmax],'k-','linewidth',2)
            
            plot(h(dimind2,dimind1),[0 threshvc(dimind2)],[threshvc(dimind2) threshvc(dimind2)],'k-','linewidth',2)
            plot(h(dimind2,dimind1),[threshvc(dimind1) threshvc(dimind1)],[0 threshvc(dimind1)],'k-','linewidth',2)
            
        else
            
            matchinds = find(strcmp(clipstrct.speclabs,templatestrct.speclabs{template_selectinds(dimind1)}));
            nonmatchinds = setdiff(find(maxinds==template_selectinds(dimind1)),matchinds);
            
            N_match = hist(clipstrct.matchmat(matchinds,template_selectinds(dimind1)),xvc)';
            N_nonmatch = hist(clipstrct.matchmat(nonmatchinds,template_selectinds(dimind1)),xvc)';
            
%             N_match = N_match / sum(N_match);
%             N_nonmatch = N_nonmatch / sum(N_nonmatch);
            
            plot(h(dimind2,dimind1),xvc,N_match,'k-o','linewidth',2)
            hold on
            plot(h(dimind2,dimind1),xvc,N_nonmatch,'k--x','linewidth',2)
            thresh = templatestrct.threshvc(template_selectinds(dimind1));
            
%             ylim([0 1])
            ylims = get(gca,'ylim');
            
            plot(h(dimind2,dimind1),[thresh-.05;thresh-.05],ylims,'k-','linewidth',2)
      
            if dimind1==1
                legend({'match','nonmatch'})
            end
            
        end
        
        xlim(xlims);
        
        if dimind2 < dimnm
            set(h(dimind2,dimind1),'xtick',[])
        end
        
        if dimind1 > 1
            set(h(dimind2,dimind1),'ytick',[])
        end
        
        if dimind2 == dimnm
            xlabel(labarr{dimind1+1});
        end
        
        if dimind1 == 1
            ylabel(labarr{dimind2+1})
        end
        
        if dimind2 == 1
            title(['error = ' num2str(round(100*templatestrct.errvc(template_selectinds(dimind1)))) '%'])
        end
        
        
        hold off
        
    end
    
    if dimnm > 1
        yticks = get(h(end,1),'xtick');
        fct = diff(get(h(1,1),'ylim')) / diff(get(h(end,1),'ylim'));
        set(h(1,1),'ytick',yticks*fct,'yticklabel',yticks)
    end
    
    
end

if dimnm == 1
   ylabel('count')
%    set(gca,'ytick',[0:.1:1],'yticklabel',[0:.1:1])
    
end

deselect_all_Callback(handles.template_axis);