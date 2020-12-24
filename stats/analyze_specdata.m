clear
clc
close all

LCbirdarr = {'RD2','SI26','PU31','BL_16'};

NEbirdarr = {'BR_2','OR_13','BR26','OR2','WH57','Y437','WH27','BRFINAL'};


% NEbirdarr = {'BRFINAL'};

% NEbirdarr = {'BR26'};


stdsal = [];
mnsal = [];
stdexp = [];
mnexp = [];
birdall = [];
sylall = [];

birdiff = zeros(numel(NEbirdarr),1);
Nexp = birdiff;
Nsal = birdiff;

for birdind = 1:numel(NEbirdarr)
    load(['./' NEbirdarr{birdind} '/specdata.mat'])
    
    %     if strcmp(NEbirdarr{birdind},'BR26')
    %         stdsal = [stdsal;specstrct_U_sal.std(1:2)'];
    %         mnsal = [mnsal;specstrct_U_sal.mn(1:2)'];
    %
    %         stdexp = [stdexp;specstrct_U_NE.std(1:2)'];
    %         mnexp = [mnexp;specstrct_U_NE.mn(1:2)'];
    %
    %
    %         sylall = [sylall;specstrct_U_sal.seq(1:2)'];
    %         birdall = [birdall;repmat(NEbirdarr(birdind),2,1)];
    %
    %     else
    %
    %         stdsal = [stdsal;specstrct_U_sal.std'];
    %         mnsal = [mnsal;specstrct_U_sal.mn'];
    %
    %         stdexp = [stdexp;specstrct_U_NE.std'];
    %         mnexp = [mnexp;specstrct_U_NE.mn'];
    %
    %
    %         sylall = [sylall;specstrct_U_sal.seq'];
    %         birdall = [birdall;repmat(NEbirdarr(birdind),numel(specstrct_U_NE.std),1)];
    %
    %     end
    
    if strcmp(NEbirdarr{birdind},'BR26')
       specstrct_U_NE.spalgnarr = specstrct_U_NE.spalgnarr(1:2); 
       specstrct_U_sal.spalgnarr = specstrct_U_sal.spalgnarr(1:2); 
       specstrct_U_sal.mn = specstrct_U_sal.mn(1:2);
       specstrct_U_NE.mn = specstrct_U_NE.mn(1:2);
       specstrct_U_sal.seq = specstrct_U_sal.seq(1:2);
    end
%     
%     if size(specstrct_U_NE.spalgnarr{1},1)<size(specstrct_U_sal.spalgnarr{1},1)
%         Ntarg = size(specstrct_U_NE.spalgnarr{1},1);
%         NU = size(specstrct_U_sal.spalgnarr{1},1);
%         
%         indsrand = randperm(NU);
%         indsrand = indsrand(1:Ntarg);
%         
%         for sylind = 1:numel(specstrct_U_sal.spalgnarr)
%             specstrct_U_sal.spalgnarr{sylind} = specstrct_U_sal.spalgnarr{sylind}(indsrand,:,:);
%         end
%         
%     end
%     
    
    [stdtmp_exp,distmp_exp,stdtot_exp,spec_exp,specdev_exp,timeinds_exp] = spec_regression(specstrct_U_NE.spalgnarr);
    [stdtmp_sal,distmp_sal,stdtot_sal,spec_sal,specdev_sal,timeinds_sal] = spec_regression(specstrct_U_sal.spalgnarr);
    
    birdiff(birdind) = median(stdtmp_sal'-stdtmp_exp');
    
    stdsal = [stdsal;stdtmp_sal'];
    mnsal = [mnsal;specstrct_U_sal.mn'];
    
    stdexp = [stdexp;stdtmp_exp'];
    mnexp = [mnexp;specstrct_U_NE.mn'];
    
    sylall = [sylall;specstrct_U_sal.seq'];
    birdall = [birdall;repmat(NEbirdarr(birdind),numel(specstrct_U_sal.seq),1)];
    
    
    Nexp(birdind) = size(specstrct_U_NE.distmat,1);
    Nsal(birdind) = size(specstrct_U_sal.distmat,1);
    
    clear specstrct_U_NE specstrct_U_sal
    
end

figure
%subplot(2,1,1)
plot(stdsal,stdexp,'o')
hold on
% plot([.6 .9],[.6 .9],'k--')
plot([.5 .8],[.5 .8],'k--')

% subplot(2,1,2)
% plot(mnsal,mnexp,'o')
% hold on
% plot([1 4],[1 4],'k--')

paramat = [stdsal stdexp mnsal mnexp];