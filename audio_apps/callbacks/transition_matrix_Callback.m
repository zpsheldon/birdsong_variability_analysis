function transition_matrix_Callback(hObject, eventdata, handles)
% hObject    handle to transition_matrix (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

 
matchstrct = get(handles.load_match_data,'UserData');
cliplabs = matchstrct.cliplabs;
seqstrct = get(handles.compile_sequences,'UserData');
 
trans_strct = get(handles.transition_matrix,'UserData');
 
minP = .000001;

labsu = seqstrct.labsu;
labsu{end} = '*';
silind = length(labsu);

 
entvc = zeros(1,4);
entvc_cont = entvc;
entvc_shuff = entvc;
entvc_cont_shuff = entvc;

 
klvc = entvc;
klvc_cont = entvc;
klvc_shuff = entvc;

 
% vcshuff = seqstrct.vc(randperm(length(seqstrct.vc)));

 
pltorder = 1;

 
for order = length(entvc):-1:1    

    
    [transmat,P,rowinds,colinds] = vc2transmatOLD(seqstrct.vc,order);

    
    transarr{order} = transmat;
    Parr{order} = P;
    rowarr{order} = rowinds;
    colarr{order} = colinds;
 
    colabs = labsu(colinds);
    
    [dmy,cont_rowinds] = find(rowinds==silind);
    cont_rowinds = setdiff(1:size(rowinds,2),cont_rowinds);
  
    cont_colinds = find(colinds~=silind);
    colabs_cont = colabs(cont_colinds);

    
    transcont = transmat(cont_rowinds,cont_colinds);
    Pcont = transcont ./ repmat(sum(transcont,2),1,size(transcont,2));
   
    if order == 1
        P_plt = P;
        Pcont_plt = Pcont;
        rowlabs_plt = seqstrct.labsu;
        colabs_plt = seqstrct.labsu;
        rowlabs_cont_plt = seqstrct.labsu(setdiff(1:length(seqstrct.labsu),silind));
        colabs_cont_plt = rowlabs_cont_plt;
        transmat_plt = transmat;
    end

    
end

 
silind = find(colarr{1}==silind);
begdist = transarr{1}(silind,:);
begP = max(begdist / sum(begdist),minP);

 
Q = sum(transarr{1}) / sum(transarr{1}(:));

 
beg_kl = sum(max(begP,minP) .* log2(max(begP ./ Q,minP)));

 
transmat = vc2transmat(seqstrct.vc(randperm(length(seqstrct.vc))),1);
begP_shuff = transmat(silind,:) / sum(transmat(silind,:));

 
beg_kl_shuff = sum(max(begP_shuff,minP) .* log2(max(begP_shuff ./ Q,minP)));

 

 
endist = transarr{1}(:,silind);
endP = max(endist / sum(endist),minP);

 
endP_shuff = transmat(:,silind)' / sum(transmat(:,silind));

 

 
end_kl = sum(max(endP',minP) .* log2(max(endP' ./ Q,minP)));
end_kl_shuff = sum(max(endP_shuff,minP) .* log2(max(endP_shuff ./ Q,minP)));

 
end_kl_nrm = end_kl - end_kl_shuff;
beg_kl_nrm = beg_kl - beg_kl_shuff;

 

 
klvc = seqkl(seqstrct.vc,4,silind);
klvc_shuff = seqkl(seqstrct.vc(randperm(length(seqstrct.vc))),4,silind);
klvc_nrm = klvc - klvc_shuff;

 
trans_strct.transarr = transarr;
trans_strct.Parr = Parr;
trans_strct.rowarr = rowarr;
trans_strct.colarr = colarr;

 
trans_strct.kl = klvc;
trans_strct.kl_nrm = klvc_nrm;
trans_strct.beg_kl = beg_kl;
trans_strct.end_kl = end_kl;
trans_strct.beg_kl_nrm = beg_kl_nrm;
trans_strct.end_kl_nrm = end_kl_nrm;

 
set(handles.transition_matrix,'UserData',trans_strct)

 
M = sum(transmat_plt,2);

 
htmp = subplot(1,1,1,'Parent',handles.analysis_panel);
cla(htmp,'reset')

 
h(1) = subplot(4,6,[2:5 [2:5]+6],'Parent',handles.analysis_panel);
bar(h(1),1:size(M,1),M(:,1));
set(h(1),'xtick',1:size(M,1),'xticklabel',rowlabs_plt','box','off');
title('syllable freq.')

 
klvc = round(100*klvc)/100;
klvc_nrm = round(100*klvc_nrm)/100;
beg_kl = round(100*beg_kl)/100;
end_kl = round(100*end_kl)/100;
beg_kl_nrm = round(100*beg_kl_nrm)/100;
end_kl_nrm = round(100*end_kl_nrm)/100;

 
entab = sprintf('%-20s%-20s%-20s\n',...
    '','kl div','nrm kl div',...
    'beg.',num2str(beg_kl),num2str(beg_kl_nrm),...
    'end.',num2str(end_kl),num2str(end_kl_nrm),...
    'seq ord. 1',num2str(klvc(1)),num2str(klvc_nrm(1)),...
    'seq ord. 2',num2str(klvc(2)),num2str(klvc_nrm(2)),...
    'seq ord. 3',num2str(klvc(3)),num2str(klvc_nrm(3)),...
    'seq ord. 4',num2str(klvc(4)),num2str(klvc_nrm(4)));

 

 
h(2) = subplot(4,6,[[2:5]+12 [2:5]+18],'Parent',handles.analysis_panel);

 
% text(0,0,entab,'horizontalalignment','center','verticalalignment','middle');
text(0,0,entab,'Interpreter','none','FontUnits','normalized');
ylim([-.7 .7])
xlim([-.5 2])

 
set(h(2),'xtick',[],'ytick',[],'box','off');

 
htmp = subplot(1,1,1,'Parent',handles.analysis_panel2);
cla(htmp,'reset')

 
h(4) = subplot(2,6,2:5,'Parent',handles.analysis_panel2);
imagesc(P_plt,[0 .999])
colorbar('peer',h(4));
set(h(4),'ydir','normal','xtick',1:size(P_plt,2),'ytick',1:size(P_plt,1),'xticklabel',colabs_plt,'yticklabel',rowlabs_plt','ylim',[.5 size(P_plt,1)+.5],'xlim',[.5 size(P_plt,2)+.5],'box','off');

 

 
h(5) = subplot(2,6,[2:5]+6,'Parent',handles.analysis_panel2);
imagesc(Pcont_plt,[0 .999])
colorbar('peer',h(5));
set(h(5),'ydir','normal','xtick',1:size(Pcont_plt,2),'ytick',1:size(Pcont_plt,1),...
    'xticklabel',colabs_cont_plt,'yticklabel',rowlabs_cont_plt','ylim',[.5 size(Pcont_plt,1)+.5],'xlim',[.5 size(Pcont_plt,2)+.5],'box','off');

 

 
% 
% trans_strct.transmat = transmat;
% trans_strct.P = P;
% trans_strct.P_nosil = Pcont;

 
trans_strct.transarr = transarr;
trans_strct.Parr = Parr;
trans_strct.rowarr = rowarr;
trans_strct.colarr = colarr;

 

 
if length(seqstrct.vc) < 10000 && ~isfield(trans_strct,'warnflg')
   warndlg('Under 10,000 clips, KL estimates above order = 2 may have large error')
   trans_strct.warnflg = 1;
end

 
set(handles.transition_matrix,'UserData',trans_strct)