clear
clc
close all

%OR2: 6/14 should be last day of song used, 734668

bird = {'BR_2','OR_13','BR26','OR2','WH57','Y437','WH27','BRFINAL','OR_46','RD_2','SI_026','PU_31','BL_16'};

bird =  {'BL_16','BR_2','BR26','BRFINAL','OR_13','OR_46','OR2','PU_31','RD_2','SI_26','WH27','WH57','Y437'};
targsyl = {'B','G','B','C','G','A','G','A','B','A','A','A','C'};
targseq = {'DB','GD','CB','CD','GF','AB','GI','AC','BG','AC','AC','AC','GC'};

birdarrDIR = {'BR_2','OR_13','BR26','OR2','WH27','OR_46'};
birdarrLC = {'RD_2','SI_026','PU_31','BL_16'};
birdarrPHE = {'OR_46','BR26','WH27','OR_2','WH57'};
birdarrNE =  {'BR_2','OR_13','BR26','OR2','WH57','Y437','WH27','BRFINAL'};

seqNdir = zeros(numel(birdarrDIR),1)+nan;
seqNundir = seqNdir;
Tundir = seqNdir;
Tdir = Tundir;

for birdind = 1:numel(birdarrDIR)
    load(['./' birdarrDIR{birdind} '/seqdata.mat'])
    
    birdind2 = find(strcmp(bird,birdarrDIR{birdind}));
    birdsyl = targsyl(birdind2);
    
    seqstrct = seqstrct_U_sal;

    sylind = find(strcmp(seqstrct.labsu,birdsyl{1}));
    
    wavdirorg = seqstrct.wavdir;
    load([wavdirorg '/wavdirinfo.mat'])
    
    [wavtms,indsort] = sort(dirstrct.daytms,'ascend');
    
    recdur = getRecDur(wavtms,[-inf inf],maxgap);
    
%     vcrl = seqstrct.vc(seqstrct.wavinds>0);
%     wavindsrl = seqstrct.wavinds(seqstrct.wavinds>0);
%     
%     cliptms = zeros(numel(vcrl),1);
%         
%     for clipind = 1:numel(vcrl)
%         cliptms(clipind) = dirstrct.daytms(dirstrct.wavinds==wavindsrl(clipind));
%     end
%     
%     [cliptms,indsort] = sort(cliptms);
%     
    
    
%     vcrl = vcrl(indsort);
    
end

sngRateUndir = seqNundir./Tundir;
sngRateDir = seqNdir./Tdir;
