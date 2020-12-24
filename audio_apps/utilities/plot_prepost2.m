function [h,lenmat,errmat,cvmat,cLmat,cUmat] = plot_prepost2(seqarr,lenarr)

h = zeros(1,4);

clrarr = {'k','b','r','m','g'};

lenmat = zeros(length(seqarr),length(lenarr));
errmat = lenmat;
cvmat = lenmat;
cvLmat = lenmat;
cvUmat = lenmat;

seqlen = zeros(1,length(lenarr));
seqlenerr = seqlen;
seqcv = seqlen;
seqcvL = seqlen;
seqcvU = seqlen;

M = length(seqarr);

for i = 1:length(lenarr)
    
    N = size(lenarr{i},1);
    
    lens = lenarr{i} ./ repmat(mean(lenarr{1}),N,1);
    
    lenmat(:,i) = mean(lens)';
    errmat(:,i) = sem(lens)';
    
    seqlens = sum(lenarr{i}')' / mean(sum(lenarr{1}')');
    seqlen(i) = mean(seqlens);
    seqlenerr(i) = sem(seqlens);
    
    
    CL = sqrt((N-1)*var(lenarr{i})/chi2inv(.025,(N-1))) ./ mean(lenarr{i});
    CU = sqrt((N-1)*var(lenarr{i})/chi2inv(.975,(N-1))) ./ mean(lenarr{i});
    
    cvmat(:,i) = cv(lenarr{i})';
    cLmat(:,i) = CL';
    cUmat(:,i) = CU';
    
    
    CL = sqrt((N-1)*var(seqlens)/chi2inv(.025,(N-1))) ./ mean(seqlens);
    CU = sqrt((N-1)*var(seqlens)/chi2inv(.975,(N-1))) ./ mean(seqlens);
    
    seqcv(i) = cv(seqlens);
    seqcvL(i) = CL;
    seqcvU(i) = CU;
    
end



figure
h(1) = subplot(2,1,1);
errorbar(lenmat,errmat)
set(gca,'xtick',1:length(seqarr),'xticklabel',seqarr,'box','off')
ylabel('norm len')

for i = 1:length(lenarr)
    lgndtxt{i} = num2str(i);
end

legend(lgndtxt);


h(2) = subplot(2,1,2);
errorbar(repmat(1:length(seqarr),length(lenarr),1)',cvmat,cvmat-cLmat,cUmat-cvmat)
set(gca,'xtick',1:length(seqarr),'xticklabel',seqarr,'box','off')
ylabel('length CV')


figure
h(3) = axes;
errorbar(lenmat',errmat')
hold on
errorbar(seqlen,seqlenerr,'k','linewidth',3)
plot([1 length(lenarr)],[1 1],'k--','linewidth',2)
ylim([.85 1.15])

figure
h(4) = axes;
errorbar(repmat(1:length(lenarr),length(seqarr),1)',cvmat',cvmat'-cLmat',cUmat'-cvmat')
hold on
errorbar(1:length(seqcv),seqcv,seqcv-seqcvL,seqcvU-seqcv,'k','linewidth',3)
