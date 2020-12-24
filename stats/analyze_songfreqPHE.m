clear
clc
close all

%OR2: 6/14 should be last day of song used, 734668

bird = {'BR_2','OR_13','BR26','OR2','WH57','Y437','WH27','BRFINAL','OR_46','RD_2','SI_026','PU_31','BL_16'};

bird =  {'BL_16','BR_2','BR26','BRFINAL','OR_13','OR_46','OR2','PU_31','RD_2','SI_26','WH27','WH57','Y437'};
targsyls = {'B','G','B','C','G','A','G','A','B','A','A','A','C'};
targseq = {'DB','GD','CB','CD','GF','AB','GI','AC','BG','AC','AC','AC','GC'};

birdarrDIR = {'BR_2','OR_13','BR26','OR2','WH27','OR_46'};
birdarrLC = {'RD_2','SI_026','PU_31','BL_16'};
birdarrPHE = {'OR_46','BR26','WH27','WH57'};
birdarrNE =  {'BR_2','OR_13','BR26','OR2','WH57','Y437','WH27','BRFINAL'};

seqNPHE = zeros(numel(birdarrPHE),1)+nan;
seqNsal = seqNPHE;
TPHE = seqNPHE;
Tsal = seqNPHE;

cutoff = 10/(24*60);

for birdind = 1:numel(birdarrPHE)
    load(['./' birdarrPHE{birdind} '/seqdata.mat'])
    
    switch birdarrPHE{birdind}
        case 'OR2'
            daylims = [-inf 734669];
        case 'WH57'
            daylims = [-inf 734656];
        otherwise
            daylims = [-inf inf];
    end
    
    birdind2 = find(strcmp(bird,birdarrPHE{birdind}));
    birdseq = targseq(birdind2);
    
    wavdirorg = seqstrct_F_sal.wavdir;
    load([wavdirorg '/wavdirinfo.mat'])
    [Tsal(birdind),wavindsval] = getRecDur(dirstrct.daytms,daylims,2/24,cutoff);
    
    vcval = seqstrct_F_sal.vc(ismember(seqstrct_F_sal.wavinds,wavindsval));
    
    birdseqinds = [find(strcmp(seqstrct_F_sal.labsu,birdseq{1}(1))) find(strcmp(seqstrct_F_sal.labsu,birdseq{1}(2)))];
    
    seqinds = findSeq(vcval,birdseqinds);
    seqNsal(birdind) = length(seqinds);
    
    
    wavdirorg = seqstrct_F_PHE.wavdir;
    load([wavdirorg '/wavdirinfo.mat'])
    [TPHE(birdind),wavindsval] = getRecDur(dirstrct.daytms,daylims,2/24,cutoff);
    
    vcval = seqstrct_F_PHE.vc(ismember(seqstrct_F_PHE.wavinds,wavindsval));
    
    birdseqinds = [find(strcmp(seqstrct_F_PHE.labsu,birdseq{1}(1))) find(strcmp(seqstrct_F_PHE.labsu,birdseq{1}(2)))];
    
    seqinds = findSeq(vcval,birdseqinds);
    seqNPHE(birdind) = length(seqinds);
    
end

Tsal = Tsal * 24 * 60;
TPHE = TPHE * 24 * 60;

sngRateSal = seqNsal./Tsal;
sngRatePHE = seqNPHE./TPHE;

figure; 

plotbird(sngRateSal,sngRatePHE,birdarrPHE,gca)
hold on
plot([0 5],[0 5],'k--')

set(gca,'box','off','ticklength',[0 0],'fontsize',12)
xlabel('songs/minute with saline','fontsize',16)
ylabel('songs/minute with PHE infus.','fontsize',16)

set(gcf,'Position',[440   570   271   228])


median(sngRatePHE-sngRateSal)

[d,p,chistat,chithresh] = calculateRateStat(seqNPHE,TPHE,seqNsal,Tsal,1000);
dstat = mean(d);
dsem = std(d);