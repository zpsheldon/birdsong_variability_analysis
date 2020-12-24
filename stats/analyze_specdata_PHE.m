clear
clc
close all

%PHEbirdarr = {'OR_46','BR26','WH57','BR_2'};

%PHEbirdarr = {'OR_46','BR26','WH57','BR_2','WH27'};
PHEbirdarr = {'OR_46','BR26','WH27','WH57'};

% PHEbirdarr = {'OR_46','BR26','BR_2','WH27'};

% 
% PHEbirdarr = {'OR_46','BR26'};

%NEbirdarr = {'BR_2','OR_13','BR26','OR2','WH57','BRFINAL'};


% NEbirdarr = {'BR26'};


stdsal = [];
mnsal = [];
stdexp = [];
mnexp = [];
birdall = [];
sylall = [];
lenall = [];

birdiff = zeros(numel(PHEbirdarr),1);
Nexp = birdiff;
Nsal = birdiff;

for birdind = 1:numel(PHEbirdarr)
    load(['./' PHEbirdarr{birdind} '/specdata_honed.mat'])
    
    [stdtmp_exp,distmp_exp,stdtot_exp,spec_exp,specdev_exp,timeinds_exp] = spec_regression(specstrct_F_PHE.spalgnarr);
    [stdtmp_sal,distmp_sal,stdtot_sal,spec_sal,specdev_sal,timeinds_sal] = spec_regression(specstrct_F_sal.spalgnarr);
    
    stdsal = [stdsal;stdtmp_sal'];
    mnsal = [mnsal;specstrct_F_sal.mn'];
    
    stdexp = [stdexp;stdtmp_exp'];
    mnexp = [mnexp;specstrct_F_PHE.mn'];
    
    for sylind = 1:numel(specstrct_F_PHE.spalgnarr)
       lenall = [lenall;size(specstrct_F_PHE.spalgnarr{sylind},3)];
    end
    
    sylall = [sylall;specstrct_F_sal.seq'];
    birdall = [birdall;repmat(PHEbirdarr(birdind),numel(specstrct_F_sal.std),1)];
    
    birdiff(birdind) = stdtot_exp-stdtot_sal;
    
    Nexp(birdind) = size(specstrct_F_PHE.distmat,1);
    Nsal(birdind) = size(specstrct_F_sal.distmat,1);

end

figure
% plot(stdsal./mnsal,stdexp./mnexp,'o')
% hold on
% plot([.1 .5],[.1 .5],'k--')

plotbird(stdsal,stdexp,birdall)
hold on
plot([.5 .8],[.5 .8],'k--')

paramat = [stdsal stdexp mnsal mnexp];