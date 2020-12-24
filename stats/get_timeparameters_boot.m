function paramat = get_timeparameters_boot(birdarr,cond)

birdiff2 = zeros(numel(birdarr),1);

stdsal = [];
mnsal = [];
seqmnsal = zeros(numel(birdarr),1);
seqstdsal = seqmnsal;
seqmnsal2 = seqmnsal;
seqmnexp2 = seqmnsal;

stdexp = [];
mnexp = [];
seqmnexp = seqmnsal;
seqstdexp = seqmnsal;

Wsal = [];
Wexp = [];
psisal = [];
psiexp = [];

Nsal = seqmnsal;
Nexp = seqmnsal;

mnsal2 = [];
mnexp2 = [];

sgall = [];

birdiff = zeros(numel(birdarr),1);

birdall = [];
birdall2 = [];

sylall = [];
sylall2 = [];

lenp = birdiff2;

for birdind = 1:numel(birdarr)
    
    birdind
 
    load(['./' birdarr{birdind} '/timedata.mat'])
    
    switch cond
        case 1
            lenseqexp = timestrct_U_NE.lenseq;
            lenseqsal = timestrct_U_sal.lenseq;
            seq = timestrct_U_sal.seqarr;
            sgall = [sgall;timestrct_U_sal.sg(:)];
            
        case 0
            lenseqexp = timestrct_U_stim.lenseq;
            lenseqsal = timestrct_U_nostim.lenseq;
            seq = timestrct_U_nostim.seqarr;
            sgall = [sgall;timestrct_U_nostim.sg(:)];
            
        case 2
            lenseqexp = timestrct_F_sal.lenseq;
            lenseqsal = timestrct_U_sal.lenseq;
            seq = timestrct_U_sal.seqarr;
            sgall = [sgall;timestrct_U_sal.sg(:)];
            
        case 3
            lenseqsal = timestrct_F_sal.lenseq;
            lenseqexp = timestrct_F_PHE.lenseq;
            seq = timestrct_F_sal.seqarr;
            sgall = [sgall;timestrct_F_sal.sg(:)];
            
    end
    
    Nsal = size(lenseqsal,1);
    Nexp = size(lenseqexp,1);

    salinds = ceil(Nsal*rand(Nsal,1));
    expinds = ceil(Nexp*rand(Nexp,1));
    
    lenseqsal = lenseqsal(salinds,:);
    lenseqexp = lenseqexp(expinds,:);
    
    lenseqsal = cleanInts(lenseqsal,5,1);
    lenseqexp = cleanInts(lenseqexp,5,1);
    
    
    birdiff2(birdind) = median(sum(lenseqsal')')-median(sum(lenseqexp')');
    lenp(birdind) = ranksum(sum(lenseqsal')',sum(lenseqexp')');
    
    
%     lenseqsal = lenseqsal(:,5:end);
%     lenseqexp = lenseqexp(:,5:end);
%     
%     
%     Q = mkQmat(seq);
%     
    Q = zeros(numel(seq),0);
    %     Q = zeros(3,0);
    
%     [Wsaltmp,psisaltmp,phisal,sigmasaltmp,omegasal] = CFAfull_spc_waitbar(lenseqsal,Q,100,2);
%     
%     [Wexptmp,psiexptmp,phisal,sigmaexptmp,omegasal] = CFAfull_spc_waitbar(lenseqexp,Q,100,2);
%     
    [Wsaltmp,psisaltmp,phisaltmp,sigmasaltmp,omegasal] = CFAfull_prior_spc(lenseqsal,Q,500,1);
    
    [Wexptmp,psiexptmp,phiexptmp,sigmaexptmp,omegasal] = CFAfull_prior_spc(lenseqexp,Q,500,1);
    
    Wsaltmp = rotate_maxaccum(Wsaltmp);
    Wexptmp = rotate_maxaccum(Wexptmp);
    
    Wsaltmp = Wsaltmp(:,1);
    Wexptmp = Wexptmp(:,1);
    
    if strcmp(birdarr{birdind},'BR26')
        Wsaltmp = Wsaltmp(end-2:end);
        Wexptmp = Wexptmp(end-2:end);
        psisaltmp = psisaltmp(end-2:end);
        psiexptmp = psiexptmp(end-2:end);
        
        lenseqexp = lenseqexp(:,end-2:end);
        lenseqsal = lenseqsal(:,end-2:end);
        
        seq = seq(end-2:end);
    end
    
%     
%     if strcmp(birdarr{birdind},'BRFINAL')
%         Wsaltmp = Wsaltmp(3:end);
%         Wexptmp = Wexptmp(3:end);
%         psisaltmp = psisaltmp(3:end);
%         psiexptmp = psiexptmp(3:end);
%         
%         lenseqexp = lenseqexp(:,3:end);
%         lenseqsal = lenseqsal(:,3:end);
%         seq = seq(3:end);
%     end

    sylall = [sylall;seq'];
    sylall2 = [sylall2;seq(2:end-1)'];
    
    Nsal(birdind) = size(lenseqsal,1);
    Nexp(birdind) = size(lenseqexp,1);
    
    if Nsal(birdind)>0 & Nexp(birdind)>0
        
        Wsal = [Wsal;Wsaltmp];
        Wexp = [Wexp;Wexptmp];
        
        psisal = [psisal;psisaltmp];
        psiexp = [psiexp;psiexptmp];
        
    end
     
    birdiff(birdind) = median(psisaltmp-psiexptmp);
    
    birdall = [birdall;repmat(birdarr(birdind),numel(Wsaltmp),1)];
%     birdall2 = [birdall2;repmat(birdarr(birdind),numel(psisaltmp)-2,1)];
    
    mnsal = [mnsal;mean(lenseqsal)'];   
    mnexp = [mnexp;mean(lenseqexp)'];

    
end


paramat = [psisal psiexp mnsal mnexp Wsal Wexp];