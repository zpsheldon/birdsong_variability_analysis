clear
clc
close all

load boot_timeparams_dirundir.mat

psidiff = zeros(bootN,1);
Wdiff = psidiff;
mndiff = psidiff;
Wper = Wdiff;
psiper = psidiff;

for bootind = 1:bootN
    
    paramat = paramarr{bootind};
    
    psidiff(bootind) = mean(paramat(:,1)-paramat(:,2));
    Wdiff(bootind) = mean(paramat(:,5)-paramat(:,6));
    mndiff(bootind) = mean(paramat(:,3)-paramat(:,4));
    
    Wper(bootind) = 100*mean((paramat(:,5)-paramat(:,6))./paramat(:,6));
    
end