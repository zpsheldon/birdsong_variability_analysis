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
birdarrPHE = {'OR_46','BR26','WH27','OR_2','WH57'};
birdarrNE =  {'BR_2','OR_13','BR26','OR2','WH57','Y437','WH27','BRFINAL'};

birdarr = birdarrLC;

durs = zeros(numel(birdarr),1);

for birdind = 1:numel(birdarr)
    load(['./' birdarr{birdind} '/seqdata.mat'])
 
    wavdirorg = seqstrct_U_nostim.wavdir;
    load([wavdirorg '/wavdirinfo.mat'])
    
    switch birdarr{birdind}
        case 'OR2'
            daylims = [-inf 734669];
        case 'WH57'
            daylims = [-inf  734656];
        case 'WH27'
            daylims = [-inf  736405];
        otherwise
            daylims = [-inf inf];
    end
    
    clear seqstrct_U_nostim
 
    dayset = dirstrct.dayset(dirstrct.dayset<dirstrct.dayset(1)+200);
    
    if(exist('seqstrct_U_stim'))
        wavdirorg = seqstrct_U_stim.wavdir;
        load([wavdirorg '/wavdirinfo.mat'])


        dayset_F = dirstrct.dayset(dirstrct.dayset<dirstrct.dayset(1)+200);


        dayset = [dayset, dayset_F]

        clear seqstrct_U_stim
    end
    
    durs(birdind) = max(dayset)-min(dayset)+1;
    
    
end
