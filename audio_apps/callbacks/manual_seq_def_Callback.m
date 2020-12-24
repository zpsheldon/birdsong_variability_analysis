function  manual_seq_def_Callback(hObject, eventdata, handles)
% hObject    handle to compile_sequences (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


min_seq_nm = str2double(get(handles.min_seq_nm,'String'));
min_seq_len = str2double(get(handles.min_seq_len,'String'));

seqstrct = get(handles.compile_sequences,'UserData');

maxrepeat = str2num(get(handles.max_repeats,'String'));
silind = find(strcmp(seqstrct.labsu,'sil'));

seq = inputdlg('Sequence to compile');
seq = seq{1};

labtmp = seqstrct.labsu;
labtmp{end} = '*';
commainds = findstr(seq,',');

if ~isempty(commainds)
    vc0 = zeros(1,length(commainds)+1);
    vc0(1) = find(strcmp(seq(1:commainds(1)-1),labtmp));
    for i = 2:length(commainds)
        vc0(i) = find(strcmp(seq(commainds(i-1)+1:commainds(i)-1),labtmp));
    end
    vc0(end) = find(strcmp(seq(commainds(end)+1:end),labtmp));
else
    vc0 = zeros(1,length(seq));
    for i = 1:length(seq)
        vc0(i) = find(strcmp(seq(i),labtmp));
    end
    
end

% [seqmat,N_seq] = vc2seqs(seqstrct.vc,min_seq_nm,length(vc0),maxrepeat,silind,vc0);

inds = findSeq(seqstrct.vc,vc0);
seqmat = vc0;
N_seq = length(inds);


% [seqmat2,N_seq2] = vc2seqs(flipud(seqstrct.vc),min_seq_nm,length(vc0),maxrepeat,silind,fliplr(vc0));

% lendiff = size(seqmat,2) - size(seqmat2,2);
% 
% if lendiff > 0
%     seqmat2 = [seqmat2 zeros(size(seqmat2,1),lendiff)];
% elseif lendiff < 0
%     seqmat = [seqmat zeros(size(seqmat,1),-lendiff)];
% end

% indsval = ~ismember(seqmat2,seqmat,'rows');
% seqmat2 = seqmat2(indsval,:);
% N_seq2 = N_seq2(indsval);
% 
% seqmatall = zeros(size(seqmat,1)+size(seqmat2,1),max(size(seqmat,2),size(seqmat2,2)));
% seqmatall(1:size(seqmat,1),1:size(seqmat,2)) = seqmat;
% seqmatall(size(seqmat,1)+1:size(seqmat,1)+size(seqmat2,1),1:size(seqmat2,2)) = fliplr(seqmat2);

% seqmat = seqmatall;
% N_seq = [N_seq;N_seq2];
seqnm = length(N_seq);

if seqnm
    M_seq = zeros(seqnm,1);
    seqtxt = [];
    
    for seqind = 1:seqnm
        seqtmp = seqstrct.labsu(seqmat(seqind,find(seqmat(seqind,:))));
        seqtmp(find(strcmp(seqtmp,'sil'))) = {'*'};
        seqs{seqind} = seqtmp;
        
        M_seq(seqind) = length(seqs{seqind});
    end
    
%     cliptot = N_seq .* M_seq;
%     [dmy,indsrt] = sort(cliptot,'descend');
%     
    [dmy,indsrt] = sort(N_seq,'descend');
    
    % [dmy,indsrt] = sort(M_seq,'descend');
    
    seqs = seqs(indsrt);
    N_seq = N_seq(indsrt);
    M_seq = M_seq(indsrt);
    seqmat = seqmat(indsrt,:);
    
    warnflg = 0;
    orgnm = seqnm;
    
    if orgnm > 100
        warnflg = 1;
        seqs = seqs(1:100);
        N_seq = N_seq(1:100);
        M_seq = M_seq(1:100);
        seqmat = seqmat(1:100,:);
    end
    
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
        warndlg([num2str(orgnm) ' sequences found, top 100 taken'])
    end
    
else
    
    h = subplot(1,1,1,'Parent',handles.analysis_panel);
    cla(h,'reset')
    title('no sequences w/ specs')
    
end