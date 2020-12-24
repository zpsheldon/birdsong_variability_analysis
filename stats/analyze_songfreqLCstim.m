clear
clc
close all

%OR2: 6/14 should be last day of song used, 734668

bird = {'BR_2','OR_13','BR26','OR2','WH57','Y437','WH27','BRFINAL','OR_46','RD_2','SI_026','PU_31','BL_16'};

bird =  {'BL_16','BR_2','BR26','BRFINAL','OR_13','OR_46','OR2','PU_31','RD_2','SI_026','WH27','WH57','Y437'};
targsyls = {'B','G','B','C','G','A','G','A','B','A','A','A','C'};
targseq = {'DB','GD','CB','CD','GF','AB','GI','AC','BG','AC','AC','AC','GC'};

birdarrDIR = {'BR_2','OR_13','BR26','OR2','WH27','OR_46'};
birdarrLC = {'RD_2','SI_026','PU_31','BL_16'};
birdarrPHE = {'OR_46','BR26','WH27','OR_2','WH57'};
birdarrNE =  {'BR_2','OR_13','BR26','OR2','WH57','Y437','WH27','BRFINAL'};

seqNdir = zeros(numel(birdarrLC),1)+nan;
seqNundir = seqNdir;
Tundir = seqNdir;
Tdir = Tundir;

daylims = [-inf inf];

for birdind = 1:numel(birdarrLC)
    load(['./' birdarrLC{birdind} '/seqdata.mat'])
    
    birdind2 = find(strcmp(bird,birdarrLC{birdind}));
    birdseq = targseq(birdind2);
    
    seqN = numel(seqstrct_U_nostim.seqs);
    
    found = 0;
    seqind = 0;
    
    while found==0 && seqind<seqN
        seqind = seqind + 1;
        isfound = find(strcmp(cell2mat(seqstrct_U_nostim.seqs{seqind}),birdseq));
        found = ~isempty(isfound & numel(seqstrct_U_nostim.seqs{seqind})==length(birdseq));
    end
    if found
        seqNundir(birdind) = seqstrct_U_nostim.N_seq(seqind);
    end

    wavdirorg = seqstrct_U_nostim.wavdir;
    load([wavdirorg '/wavdirinfo.mat'])
    
    Tundir(birdind) = getRecDur(dirstrct.daytms,daylims,3/24);
    
    while found==0 && seqind<seqN
        seqind = seqind + 1;
        isfound = find(strcmp(cell2mat(seqstrct_U_stim.seqs{seqind}),birdseq));
        found = ~isempty(isfound & numel(seqstrct_U_stim.seqs{seqind})==length(birdseq));
    end
    if found
        seqNdir(birdind) = seqstrct_U_stim.N_seq(seqind);
    end

    wavdirorg = seqstrct_U_stim.wavdir;
    load([wavdirorg '/wavdirinfo.mat'])
    Tdir(birdind) = getRecDur(dirstrct.daytms,daylims,3/24);
 
end

Tundir = Tundir * 24 * 60;
Tdir = Tdir * 24 * 60;

sngRateUndir = seqNundir./Tundir;
sngRateDir = seqNdir./Tdir;


% figure; plot(sngRateUndir,sngRateDir,'o')
% hold on
% plot([0 2],[0 2],'k--')

median(sngRateDir-sngRateUndir)

[d,p,chistat,chithresh] = calculateRateStat(seqNdir,Tdir,seqNundir,Tundir,1000);
dstat = mean(d);
dsem = std(d);

figure; 
plotbird(sngRateUndir,sngRateDir,birdarrLC,gca)
legend(birdarrLC)
hold on
plot([0 3],[0 3],'k--')

set(gca,'box','off','ticklength',[0 0],'fontsize',12)
xlabel('songs/minute without LC stim.','fontsize',16)
ylabel('songs/minute with LC stim.','fontsize',16)

set(gcf,'Position',[440   570   271   228])


