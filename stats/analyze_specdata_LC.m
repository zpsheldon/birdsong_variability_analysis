clear
clc
close all

LCbirdarr = {'RD_2','SI_026','PU_31','BL_16'};
%LCbirdarr = {'RD2','SI26','PU31','BL_16'};
LCbirdarr2 = {'RD2','SI026','PU31','BL16'};

%NEbirdarr = {'BR_2','OR_13','BR26','OR2','WH57','BRFINAL'};


% NEbirdarr = {'BR26'};


stdsal = [];
mnsal = [];
stdexp = [];
mnexp = [];
birdall = [];
sylall = [];

birdiff = zeros(numel(LCbirdarr),1);
Nexp = birdiff;
Nsal = birdiff;

LCbirdarr_all = {};

for birdind = 1:numel(LCbirdarr)
    load(['./' LCbirdarr{birdind} '/specdata.mat'])
    
    [stdtmp_exp,distmp_exp] = spec_regression(specstrct_U_stim.spalgnarr);
    [stdtmp_sal,distmp_sal] = spec_regression(specstrct_U_nostim.spalgnarr);
    
    stdsal = [stdsal;stdtmp_sal'];
    mnsal = [mnsal;specstrct_U_nostim.mn'];
    
    stdexp = [stdexp;stdtmp_exp'];
    mnexp = [mnexp;specstrct_U_stim.mn'];
    
    sylall = [sylall;specstrct_U_stim.seq'];
    birdall = [birdall;repmat(LCbirdarr(birdind),numel(specstrct_U_stim.std),1)];
    
    Nexp(birdind) = size(specstrct_U_stim.distmat,1);
    Nsal(birdind) = size(specstrct_U_nostim.distmat,1);

end

figure
h = plotbird(stdsal,stdexp,birdall)
hold on
plot([.5 .9],[.5 .9],'k--')
set(gca,'ylim',[.5,.8])
set(gca,'xlim',[.5,.8])


paramat = [stdsal stdexp mnsal mnexp];