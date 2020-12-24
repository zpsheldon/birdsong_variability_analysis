clear
clc
% close all

birdarr = {'BR_2','OR_13','BR26','OR2','OR_46','WH27','WH57'};
%PHEbirdarr = {'OR_46','BR26','WH27'};

% birdarr = {'BR_2','OR_13','BR26','WH57','OR_46','WH27'};
% 
% birdarr = {'OR2'};

%NEbirdarr = {'BR_2','OR_13','BR26','OR2','WH57','BRFINAL'};


% NEbirdarr = {'BR26'};


stdsal = [];
mnsal = [];
stdexp = [];
mnexp = [];
birdall = [];
sylall = [];
lenall = [];

stdexptot = zeros(numel(birdarr),1);
stdsaltot = stdexptot;

stdexptot2 = zeros(numel(birdarr),1);
stdsaltot2 = stdexptot2;

birdiff = zeros(numel(birdarr),1);
Nexp = birdiff;
Nsal = birdiff;

freqinds = 1:99;

for birdind = 1:numel(birdarr)
    
    if exist(['./' birdarr{birdind} '/specdata_honed.mat'])
        
        load(['./' birdarr{birdind} '/specdata_honed.mat'])
        
    else
        load(['./' birdarr{birdind} '/specdata.mat'])
    end
    
    if size(specstrct_F_sal.spalgnarr{1},1)<size(specstrct_U_sal.spalgnarr{1},1)
        'here'
        Ntarg = size(specstrct_F_sal.spalgnarr{1},1);
        NU = size(specstrct_U_sal.spalgnarr{1},1);
        
        indsrand = randperm(NU);
        indsrand = indsrand(1:Ntarg);
        
        for sylind = 1:numel(specstrct_U_sal.spalgnarr)
            specstrct_U_sal.spalgnarr{sylind} = specstrct_U_sal.spalgnarr{sylind}(indsrand,:,:);
        end
        
    end
    
    
    [stdtmp_exp,distmp_exp,stdtot_exp,spec_exp,specdev_exp,timeinds_exp] = spec_regression(specstrct_F_sal.spalgnarr);
    [stdtmp_sal,distmp_sal,stdtot_sal,spec_sal,specdev_sal,timeinds_sal] = spec_regression(specstrct_U_sal.spalgnarr);
    
%     stdtmp_exp = specstrct_F_sal.std;
%     stdtmp_sal = specstrct_U_sal.std;
    
    stdsal = [stdsal;stdtmp_sal'];
    mnsal = [mnsal;specstrct_U_sal.mn'];
    
    stdexp = [stdexp;stdtmp_exp'];
    mnexp = [mnexp;specstrct_F_sal.mn'];
    
    sylall = [sylall;specstrct_U_sal.seq'];
    birdall = [birdall;repmat(birdarr(birdind),numel(specstrct_U_sal.std),1)];
    
    Nexp(birdind) = size(specstrct_F_sal.distmat,1);
    Nsal(birdind) = size(specstrct_U_sal.distmat,1);
    
    %birdiff(birdind) = median(specstrct_U_sal.std'-specstrct_F_sal.std');
    
    birdiff(birdind) = median(stdtmp_sal-stdtmp_exp);
    
    stdexptot(birdind) = stdtot_exp;
    stdsaltot(birdind) = stdtot_sal;
    
    stdexptot2(birdind) = median(stdtmp_exp);
    stdsaltot2(birdind) = median(stdtmp_sal);
    
    for sylind = 1:numel(specstrct_F_sal.spalgnarr)
        lenall = [lenall;size(specstrct_F_sal.spalgnarr{sylind},3)];
    end
    
end

cvsal = stdsal./mnsal;
cvexp = stdexp./mnexp;

minval = min(stdsal);
maxval = max(stdsal);

figure
plot(stdsal,stdexp,'o')
hold on
plot([minval maxval],[minval maxval],'k--')

paramat = [stdsal stdexp mnsal mnexp];