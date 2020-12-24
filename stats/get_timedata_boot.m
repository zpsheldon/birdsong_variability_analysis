clear
clc
close all

bootN = 500;

% 
% birdarr = {'RD_2','SI_026','PU_31','BL_16'};
% cond = 0;

birdarr = {'BR_2','OR_13','BR26','OR2','WH27','OR_46','WH57'};
cond = 2;

% birdarr =  {'BR_2','OR_13','BR26','OR2','WH57','Y437','WH27','BRFINAL'};
% 
% cond = 1;
% % 
% birdarr = {'Y437'};

% BR26: have data and anl HONED~, 2
% BR2: have data and anl HONED*, 5
% OR46: have data and anl HONED*, 6
% WH57: have data and anl HONED*, 5
% OR2: have data and anl, HONED~, 4
% OR13: have data and anl, HONED, M 6
% birdarr = {'WH57'};

birdN = numel(birdarr);

paramarr = cell(bootN,1);
psidiff = zeros(bootN,1);
Wdiff = psidiff;
mndiff = psidiff;

for bootind = 1:bootN
    
    bootind
    
    birdinds = ceil(birdN*rand(birdN,1));
    birdboot = birdarr(birdinds);
    
    paramat = get_timeparameters_boot(birdboot,cond);
    
    paramarr{bootind} = paramat;
    
   % paramat = [psisal psiexp mnsal mnexp Wsal Wexp];
   
   psidiff(bootind) = mean(paramat(:,1)-paramat(:,2));
   Wdiff(bootind) = mean(paramat(:,5)-paramat(:,6));
   mndiff(bootind) = mean(paramat(:,3)-paramat(:,4));

end