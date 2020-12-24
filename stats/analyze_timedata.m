clear
clc
close all

birdarr = {'RD_2','SI_026','PU_31','BL_16'};

cond = 0;

% birdarr = {'OR_46','BR26','WH27','WH57'};
% cond = 3;
% 
% birdarr =  {'BR_2','OR_13','BR26','OR2','WH57','Y437','WH27','BRFINAL'};
% 
% cond = 1;

birdarr = {'BR_2','OR_13','BR26','OR2','WH27','OR_46','WH57'};
cond = 2;
% % 
% birdarr = {'Y437'};

% BR26: have data and anl HONED~, 2
% BR2: have data and anl HONED*, 5
% OR46: have data and anl HONED*, 6
% WH57: have data and anl HONED*, 5
% OR2: have data and anl, HONED~, 4
% OR13: have data and anl, HONED, M 6
% birdarr = {'WH57'};

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
    
%     if cond==2
%         load(['./' birdarr{birdind} '/timedata.mat'])
%     else
%         load(['./' birdarr{birdind} '/timedata.mat'])
%     end
    
%     if exist(['./' birdarr{birdind} '/timedata_honed.mat'])
%         load(['./' birdarr{birdind} '/timedata_honed.mat'])
%     else
%         load(['./' birdarr{birdind} '/timedata.mat'])
%     end
    
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
    
 %   birdiff(birdind) = median(psisaltmp(2:end-1)./mean(lenseqsal(:,2:end-1))'-psiexptmp(2:end-1)./mean(lenseqexp(:,2:end-1))');
    
    
    birdiff(birdind) = median(psisaltmp-psiexptmp);
    
    birdall = [birdall;repmat(birdarr(birdind),numel(Wsaltmp),1)];
%     birdall2 = [birdall2;repmat(birdarr(birdind),numel(psisaltmp)-2,1)];
    
    mnsal = [mnsal;mean(lenseqsal)'];
    stdsal = [stdsal;std(lenseqsal)'];
    
%     mnsal2 = [mnsal2;mean(lenseqsal(:,2:end-1))'];
%     mnexp2 = [mnexp2;mean(lenseqexp(:,2:end-1))'];
    
    mnexp = [mnexp;mean(lenseqexp)'];
    stdexp = [stdexp;std(lenseqexp)'];
    
    seqmnsal(birdind) = mean(sum(lenseqsal')')';
    seqmnexp(birdind) = mean(sum(lenseqexp')')';
    
%     seqmnsal2(birdind) = mean(sum(lenseqsal(:,2:end-1)')')';
%     seqmnexp2(birdind) = mean(sum(lenseqexp(:,2:end-1)')')';
    
    seqstdsal(birdind) = std(sum(lenseqsal')')';
    seqstdexp(birdind) = std(sum(lenseqexp')')';
    
    seqstdsal(birdind) = sqrt(sum(psisaltmp')')';
    seqstdexp(birdind) = sqrt(sum(psiexptmp')')';
    

   % clear timestrct_U_sal timestrct_U_NE
    
end

cvsal = seqstdsal./seqmnsal;
cvexp = seqstdexp./seqmnexp;

% figure
% plot(cvsal,cvexp,'o')
% hold on
% plot([.01 .1],[.01 .1],'k--')

figure;
subplot(2,1,1)
plot(Wsal,Wexp,'o'); hold on; plot([-1 4],[-1 4],'k--'); xlim([-1 4]); ylim([-1 4])

subplot(2,1,2)
plot(sqrt(psisal),sqrt(psiexp),'o'); hold on; plot([0 5],[0 5],'k--'); xlim([0 5]); ylim([0 5])


figure;
subplot(3,1,1)
plot(mnsal,mnexp,'o'); hold on;
plot([30 200],[30 200],'k--'); % xlim([-.005 .1]); ylim([-.005 .1])

subplot(3,1,2)
plot(Wsal./mnsal,Wexp./mnexp,'o'); hold on;
plot([-.005 .1],[-.005 .1],'k--'); % xlim([-.005 .1]); ylim([-.005 .1])

subplot(3,1,3)
plot(sqrt(psisal)./mnsal,sqrt(psiexp)./mnexp,'o'); hold on;
plot([0 .1],[0 .1],'k--');% xlim([-.005 .1]); ylim([-.005 .1])


paramat = [psisal psiexp mnsal mnexp Wsal Wexp];