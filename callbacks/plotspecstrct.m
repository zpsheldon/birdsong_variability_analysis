function plotspecstrct(handles,axis_handle,slider_handle,specstrct,xlen,ylen,srtopt)

if isempty(specstrct.specarr)
    
    cla(axis_handle)
    
else
    
    set(handles.mainfig,'CurrentAxes',axis_handle);
    propstrct = get(axis_handle,'UserData');
    
    ylims = [];
    if isfield(propstrct,'h')
        ylims = get(axis_handle,'ylim');
    end
    
    cla(axis_handle,'reset');
    
    xbins = round(xlen / (propstrct.f_winadv * 1000 / propstrct.fs));
    specspc_x = round(50 / (propstrct.f_winadv * 1000 / propstrct.fs));
    
    propstrct.labeled = zeros(1,length(specstrct.specarr));
    propstrct.ytops = propstrct.labeled;
    propstrct.xpos = propstrct.labeled;
    
    specspc_y = propstrct.f_winlen / 5;
    
    specnm = length(specstrct.specarr);
    border = round(10 / (propstrct.f_winadv * 1000 / propstrct.fs));
    
    
    x0 = border;
    y0 = 0;
    rows = 1;
    
    ampmax = 0;
    ampmin = 0;
    
    propstrct.h = zeros(specnm,1);
    propstrct.selectvc = propstrct.h;
    
    switch srtopt
        case 'label'
            
            if isfield(specstrct,'speclens')
                
                labvc = zeros(specnm,1);
                for specindind = 1:specnm
                    labvc(specindind) = specstrct.speclabs{specindind}(1)+0;
                end
                
                [dmy,srtinds] = sortrows([labvc -specstrct.speclens'],[1 2]);
                
            else
                [dmy,srtinds] = sort(specstrct.speclabs);
                
            end
            
        case 'length'
            [dmy,srtinds] = sort(specstrct.speclens,'descend');
    end
    
    hold on
    
    propstrct.speclabs = {};
    
    for specindind = 1:specnm
        
        specind = srtinds(specindind);
        
        speclen = size(specstrct.specarr{specind},2);
        
        if x0 + speclen - 1 > xbins - border
            y0 = y0 - propstrct.f_winlen - specspc_y + 1;
            rows = rows + 1;
            x0 = border;
        end
        
        xinds = x0:x0+speclen-1;
        yinds = y0-propstrct.f_winlen+1:y0;
        specmax = max(specstrct.specarr{specind}(:));
        specmin = min(specstrct.specarr{specind}(:));
        
        
        if specmax > ampmax
            ampmax = specmax;
        end
        
        if specmin < ampmin
            ampmin = specmin;
        end
        

        propstrct.h(specind) = imagesc(xinds,yinds,specstrct.specarr{specind});
        
        if isfield(specstrct,'featms') & length(specstrct.featms) >= specind & ~isempty(specstrct.featms{specind})
            featms = specstrct.featms{specind};
            xinds2 = repmat(xinds(featms),2,1);
            yinds2 = repmat([min(yinds);max(yinds)],1,length(featms));
            plot(xinds2,yinds2,'k-','linewidth',1.5)
        end
        
        propstrct.ytops(specind) = yinds(end);
        propstrct.xpos(specind) = xinds(round(length(xinds)/2));
        propstrct.speclabs{specind} = specstrct.speclabs{specind};
        
%         if isfield(specstrct,'scr')
% %             propstrct.speclabs{specind} = [specstrct.speclabs{specind} '-' num2str(round(100*(specstrct.scr(specind)))/100)];
%             propstrct.speclabs{specind} = specstrct.speclabs{specind};
%         elseif isfield(specstrct,'threshvc')
%             propstrct.speclabs{specind} = [specstrct.speclabs{specind} '-' num2str(round(100*(specstrct.threshvc(specind)))/100)];
%         else
%             propstrct.speclabs{specind} = specstrct.speclabs{specind};
%         end
        
        set(propstrct.h(specind),'ButtonDownFcn',handles.clickfun,'UserData',specind);
        
        x0 = x0 + speclen + specspc_x;
        
    end
    
    propstrct.ymax = propstrct.f_winlen/4;
    propstrct.ymin = y0;
    
    propstrct.srtopt = srtopt;
    
        
    if ylen == 0
        ylen = rows;        
    end
    propstrct.yscope = (propstrct.f_winlen + specspc_y) * ylen + specspc_y;
    
    
    ylims = [propstrct.ymax-propstrct.yscope propstrct.ymax];
    
    set(axis_handle,'clim',.9*[ampmin ampmax],'xlim',[0 xbins])
    set(axis_handle,'ylim',ylims)
    
    set(axis_handle,'UserData',propstrct,'ytick',[])
    
    if ~isempty(slider_handle)
        %         set(slider_handle,'SliderStep',[min(.5,1/rows) min(.5,1/rows)])
        clip_slider_Callback(slider_handle, 0, handles);
    else
        
        for clipind = 1:length(propstrct.speclabs)
            text_x = propstrct.xpos(clipind);
            text_y = propstrct.ytops(clipind);
            text(text_x,text_y,propstrct.speclabs(clipind),'HorizontalAlignment','center','VerticalAlignment','bottom','Interpreter','none');
        end
        
    end
    
end

