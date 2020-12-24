function plot_seq_struct_Callback(hObject, eventdata, handles)

seqstrct = get(handles.compile_sequences,'UserData');
propstrct = get(get(handles.analysis_panel,'Children'),'UserData');
matchstrct = get(handles.load_match_data,'UserData');
transtrct = get(handles.transition_matrix,'UserData');

seqind = find(propstrct.selectvc);
seq = seqstrct.seqs{seqind};
seqvc = seqstrct.seqmat(seqind,:);
seqvc = seqvc(find(seqvc));

% [Nvc,Pvc,Nvcb,Pvcb] = calcseqstats(seqvc,seqstrct.vc);

[Nvc,Pvc,Nvcb,Pvcb] = calcseqstats(seqvc,seqstrct.vc);
transmat1 = transtrct.transarr{1};
transmat2 = transtrct.transarr{2};

sylnm = length(seq);
Pvc1 = zeros(1,sylnm-1);
Pvc2 = zeros(1,sylnm-2);
% Pvcb1 = Pvc1;

Pvc1_l = Pvc1;
Pvc1_u = Pvc1;

Pvc2_l = Pvc2;
Pvc2_u = Pvc2;

% Pvcb1_l = Pvc1;
% Pvcb1_u = Pvc1;

Pvc_l = Pvc1;
Pvc_u = Pvc1;
% Pvcb_l = Pvc1;
% Pvcb_u = Pvc1;

seqstrct.labsu(find(strcmp(seqstrct.labsu,'sil'))) = {'*'};

for i = 1:sylnm-1
    
    [Pvc_l(i),Pvc_u(i)] = propCI(Pvc(i),Nvc(i),.05);
    ind1 = find(strcmp(seqstrct.labsu(transtrct.rowarr{1}),seq{i}));
    ind2 = find(strcmp(seqstrct.labsu(transtrct.colarr{1}),seq{i+1}));
    
    Pvc1(i) = transmat1(ind1,ind2) / sum(transmat1(ind1,:));
    [Pvc1_l(i),Pvc1_u(i)] = propCI(Pvc1(i),sum(transmat1(ind1,:)),.05);
    
    if i > 1
        ind1 = find(strcmp(seqstrct.labsu(transtrct.rowarr{2}(1,:)),seq{i-1}) & strcmp(seqstrct.labsu(transtrct.rowarr{2}(2,:)),seq{i}));    
        ind2 = find(strcmp(seqstrct.labsu(transtrct.colarr{2}),seq{i+1}));
        Pvc2(i-1) = transmat2(ind1,ind2) / sum(transmat2(ind1,:));
        [Pvc2_l(i-1),Pvc2_u(i-1)] = propCI(Pvc2(i-1),sum(transmat2(ind1,:)),.05);
    end
 
%     Pvcb1(i) = transmat(ind1,ind2) / sum(transmat(:,ind2));
    
%     [Pvcb1_l(i),Pvcb1_u(i)] = propCI(Pvcb1(i),sum(transmat(ind2,:)),.05);
    
%     [Pvcb_l(i),Pvcb_u(i)] = propCI(Pvcb(i),Nvcb(i),.05);
    
end

Nvc = Nvc / matchstrct.wavnm;
% Nvcb = Nvcb / matchstrct.wavnm;

h = subplot(1,1,1,'Parent',handles.analysis_panel2);
cla(h,'reset');

h = [0 0 0];
h(1) = subplot(2,1,1,'Parent',handles.analysis_panel2);
plot([1:length(Nvc)]-.5,Nvc,'b-o');
% hold on
% plot([1:length(Nvc)]-.5,Nvcb,'r-x');
ylabel('mean # per .wav file')
ylim([0 max(Nvc)*1.1])

set(h(1),'xtick',[1:length(seq)]-.5,'xticklabel',[],'xlim',[0 length(seq)],'box','off')

h(2) = subplot(2,1,2,'Parent',handles.analysis_panel2);
errorbar(1:length(Pvc),Pvc,Pvc-Pvc_l,Pvc_u-Pvc,'b-');
hold on
errorbar(1:length(Pvc1),Pvc1,Pvc1-Pvc1_l,Pvc1_u-Pvc1,'r-');
hold on
errorbar(2:length(Pvc2)+1,Pvc2,Pvc2-Pvc2_l,Pvc2_u-Pvc2,'k-');
ylim([0 1])
ylabel('forward trans. prob.')

legend({'conditional','1st order','2nd order'},'Location','best')
set(h(2),'xtick',[1:length(seq)]-.5,'xticklabel',seq,'xlim',[0 length(seq)],'box','off')

