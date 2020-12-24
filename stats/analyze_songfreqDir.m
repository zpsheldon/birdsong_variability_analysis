clear
clc
close all

%OR2: 6/14 should be last day of song used, 734668
%WH 57: last should be 734656

bird = {'BR_2','OR_13','BR26','OR2','WH57','Y437','WH27','BRFINAL','OR_46','RD_2','SI_026','PU_31','BL_16'};

bird =  {'BL_16','BR_2','BR26','BRFINAL','OR_13','OR_46','OR2','PU_31','RD_2','SI_26','WH27','WH57','Y437'};
targsyls = {'B','G','B','C','G','A','G','A','B','A','A','A','C'};
targseq = {'DB','GD','CB','CD','GF','AB','GI','AC','BG','AC','AC','AC','GC'};

birdarrDIR = {'BR_2','OR_13','BR26','OR2','WH27','OR_46','WH57'};
birdarrLC = {'RD_2','SI_026','PU_31','BL_16'};
birdarrPHE = {'OR_46','BR26','WH27','OR_2','WH57'};
birdarrNE =  {'BR_2','OR_13','BR26','OR2','WH57','Y437','WH27','BRFINAL'};

seqNdir = zeros(numel(birdarrDIR),1)+nan;
seqNundir = seqNdir;
Tundir = seqNdir;
Tdir = Tundir;

cutoff = 10/(24*60);
% cutoff = inf;

for birdind = 1:numel(birdarrDIR)
    load(['./' birdarrDIR{birdind} '/seqdata.mat'])
    
    birdind2 = find(strcmp(bird,birdarrDIR{birdind}));
    birdseq = targseq(birdind2);
    
    seqN = numel(seqstrct_U_sal.seqs);
    
    found = 0;
    seqind = 0;
    
    while found==0 && seqind<seqN
        seqind = seqind + 1;
        isfound = find(strcmp(cell2mat(seqstrct_U_sal.seqs{seqind}),birdseq));
        found = ~isempty(isfound & numel(seqstrct_U_sal.seqs{seqind})==length(birdseq));
    end
    if found
        seqNundir(birdind) = seqstrct_U_sal.N_seq(seqind);
    end
    
    wavdirorg = seqstrct_U_sal.wavdir;
    load([wavdirorg '/wavdirinfo.mat'])
    
    switch birdarrDIR{birdind}
        case 'OR2'
            daylims = [-inf 734669];
        case 'WH57'
            daylims = [-inf 734656];
        otherwise
            daylims = [-inf inf];
    end
    
    
    Tundir(birdind) = getRecDur(dirstrct.daytms,daylims,3/24);
    
%     
%     wavdirorg = seqstrct_U_sal.wavdir;
%     load([wavdirorg '/wavdirinfo.mat'])
%     
%     [Tundir(birdind),wavindsval] = getRecDur(dirstrct.daytms,daylims,2/24,cutoff);
%     
%     vcval = seqstrct_U_sal.vc(ismember(seqstrct_U_sal.wavinds,wavindsval));
%     
%     birdseqinds = [find(strcmp(seqstrct_U_sal.labsu,birdseq{1}(1))) find(strcmp(seqstrct_U_sal.labsu,birdseq{1}(2)))];
%     
%     seqinds = findSeq(vcval,birdseqinds);
%     seqNundir(birdind) = length(seqinds);
    
    
    wavdirorg = seqstrct_F_sal.wavdir;
    load([wavdirorg '/wavdirinfo.mat'])
    [Tdir(birdind),wavindsval] = getRecDur(dirstrct.daytms,daylims,3/24,cutoff);
    
    vcval = seqstrct_F_sal.vc(ismember(seqstrct_F_sal.wavinds,wavindsval));
    
    birdseqinds = [find(strcmp(seqstrct_F_sal.labsu,birdseq{1}(1))) find(strcmp(seqstrct_F_sal.labsu,birdseq{1}(2)))];
    
    seqinds = findSeq(vcval,birdseqinds);
    seqNdir(birdind) = length(seqinds);
    %
    %     while found==0 && seqind<seqN
    %         seqind = seqind + 1;
    %         isfound = find(strcmp(cell2mat(seqstrct_F_sal.seqs{seqind}),birdseq));
    %         found = ~isempty(isfound & numel(seqstrct_F_sal.seqs{seqind})==length(birdseq));
    %     end
    %     if found
    %         seqNdir(birdind) = seqstrct_F_sal.N_seq(seqind);
    %     end
    
    
end

Tundir = Tundir * 24 * 60;
Tdir = Tdir * 24 * 60;

sngRateUndir = seqNundir./Tundir;
sngRateDir = seqNdir./Tdir;

[d,p,chistat,chithresh] = calculateRateStat(seqNdir,Tdir,seqNundir,Tundir,1000);
dstat = mean(d);
dsem = std(d);


figure
plotbird(sngRateUndir,sngRateDir,birdarrDIR,gca)
legend(birdarrDIR)
hold on
plot([0 5],[0 5],'k--')

set(gca,'box','off','ticklength',[0 0],'fontsize',12)
xlabel('songs/minute without female','fontsize',16)
ylabel('songs/minute with female','fontsize',16)

set(gcf,'Position',[440   570   271   228])

% d = sngRateDir-sngRateUndir;
% 
% [d,p] = calculateRateStat(seqNdir,Tdir,seqNundir,Tundir,1000);
% dstat = median(d);
% 
% sngRateMed = median(sngRateDir-sngRateUndir);
% sngRateMedSEM = std(bootstrp(500,@median,sngRateDir-sngRateUndir));
% 
% figure
% bar(1,sngRateMed)
% hold on
% errorbar(1,sngRateMed,sngRateMedSEM);


