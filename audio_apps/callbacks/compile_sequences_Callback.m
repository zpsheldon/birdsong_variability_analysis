function compile_sequences_Callback(hObject, eventdata, handles)
% hObject    handle to compile_sequences (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


min_seq_nm = str2double(get(handles.min_seq_nm,'String'));
min_seq_len = str2double(get(handles.min_seq_len,'String'));

seqstrct = get(handles.compile_sequences,'UserData');

maxrepeat = str2num(get(handles.max_repeats,'String'));
silind = find(strcmp(seqstrct.labsu,'sil'));
[seqmat,N_seq] = vc2seqs(seqstrct.vc,min_seq_nm,min_seq_len,maxrepeat,silind);
seqnm = length(N_seq);

if seqnm
    M_seq = zeros(seqnm,1);
    Mu_seq = M_seq;
    seqtxt = [];
    
    for seqind = 1:seqnm
        seqtmp = seqstrct.labsu(seqmat(seqind,find(seqmat(seqind,:))));
        seqtmp(find(strcmp(seqtmp,'sil'))) = {'*'};
        seqs{seqind} = seqtmp;
        
        M_seq(seqind) = length(seqs{seqind});
        Mu_seq(seqind) = length(unique(seqs{seqind}));
    end
    
%     cliptot = N_seq .* Mu_seq;
%     [dmy,indsrt] = sort(cliptot,'descend');
%     
    [dmy,indsrt] = sort(N_seq,'descend');
    % [dmy,indsrt] = sort(M_seq,'descend');
    
    seqs = seqs(indsrt);
    N_seq = N_seq(indsrt);
    M_seq = M_seq(indsrt);
    seqmat = seqmat(indsrt,:);
    
    warnflg = 0;
    
    if seqnm > 100
        warnflg = 1;
        seqs = seqs(1:100);
        N_seq = N_seq(1:100);
        M_seq = M_seq(1:100);
        seqmat = seqmat(1:100,:);
    end
    
    %     [seqmat,indsrt] = sortrows(seqmat);
    %     seqs = seqs(indsrt);
    %     N_seq = N_seq(indsrt);
    %     M_seq = M_seq(indsrt);
    %
    %     cliptot = N_seq .* M_seq;
    %
    %     warnflg = 0;
    %
    %     if seqnm > 100
    %         warnflg = 1;
    %         [totsrt,indsrt] = sort(cliptot,'descend');
    %         cutoff = totsrt(100);
    %         indsval = find(cliptot>=cutoff);
    %         indsval = indsval(1:100);
    %         seqs = seqs(indsval);
    %         N_seq = N_seq(indsval);
    %         M_seq = M_seq(indsval);
    %         seqmat = seqmat(indsval,:);
    %     end
    
    digstr = num2str(floor(log10(max(N_seq)))+1);
    
    charnm = 0;
    for seqindtmp = 1:length(seqs)
        lentmp = length(cell2mat(seqs{seqindtmp}));
        if lentmp > charnm
            charnm = lentmp;
        end
    end
    charnm = num2str(charnm+1);
    
    for seqindtmp = 1:length(seqs)
        seqtxt = [seqtxt;sprintf(['%' charnm 's: %' digstr 'd'], cell2mat(seqs{seqindtmp}), N_seq(seqindtmp))];
    end
    
    colspc = 20;
    
    h = subplot(1,1,1,'Parent',handles.analysis_panel);
    cla(h,'reset')
    
    [txtrows,txtcols] = size(seqtxt);
    
    maxcols = 3;
    maxrows = ceil(txtrows / maxcols);
    
    k = 1;
    
    seqstrct.h_seq = zeros(txtrows,1);
    
    for i = 1:maxcols
        for j = 1:maxrows
            
            if k <= size(seqtxt,1)
                seqstrct.h_seq(k) = ...
                    text(i*(txtcols+colspc),(maxrows-j+1)/maxrows,seqtxt(k,:),'HorizontalAlignment','right',...
                    'VerticalAlignment','top','FontSize',10);
                
                set(seqstrct.h_seq(k),'ButtonDownFcn',handles.clickfun,'UserData',k);
                k = k + 1;
            end
            
        end
    end
    
    if txtcols > 0
        
        xlim([0 maxcols*(txtcols+colspc)])
        ylim([0 1])
        set(gca,'xtick',[],'ytick',[],'box','off')
        
        seqstrct.seqs = seqs;
        seqstrct.N_seq = N_seq;
        seqstrct.M_seq = M_seq;
        seqstrct.seqmat = seqmat;
        
        title('frequent sequences')
        
    else
        title('no sequences w/ specs')
    end
    
    set(handles.compile_sequences,'UserData',seqstrct)
    
    if warnflg
        warndlg([num2str(seqnm) ' sequences found, top 100 taken'])
    end
    
else
    
    h = subplot(1,1,1,'Parent',handles.analysis_panel);
    cla(h,'reset')
    title('no sequences w/ specs')
    
end