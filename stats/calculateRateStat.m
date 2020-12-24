function [d,p,chistat,chithresh] = calculateRateStat(seqNdir,Tdir,seqNundir,Tundir,bootN)

birdN = numel(seqNdir);
% 
% O = [seqNdir seqNundir];
% 
% E = repmat(((seqNdir + seqNundir)./(Tdir + Tundir)),1,2).*([Tdir Tundir]);
% 
% R1 = seqNdir./Tdir;
% R2 = seqNundir./Tundir;
Rnull = (seqNdir+seqNundir)./(Tdir+Tundir);

chistat = zeros(bootN,1);
chithresh = chistat;
d = chistat;

for bootind = 1:bootN
    birdinds = ceil(birdN*rand(birdN,1));
    
    seqNdirboot = seqNdir(birdinds);
    seqNundirboot = seqNundir(birdinds);
    
    Tdirboot = Tdir(birdinds);
    Tundirboot = Tundir(birdinds);
    
    Rnullboot = Rnull(birdinds);
    
    seqNdirexp = Tdirboot.*Rnullboot;
    seqNundirexp = Tundirboot.*Rnullboot;
%     
%     chistatmp = ((seqNdirboot-seqNdirexp).^2)./seqNdirexp+((seqNundirboot-seqNundirexp).^2)./seqNundirexp;
%     chistat(bootind) = sum(chistatmp);
%     chithresh(bootind) = chi2cdf(chistat(bootind),sum(seqNdirexp+seqNundirexp)-1);
%     
% %     chi2inv(.05,100);
%     
    d(bootind) = (sum(seqNdirboot)/sum(Tdirboot)./(sum(seqNundirboot)/sum(Tundirboot))-1);
    
end

p = sum(d<0)/bootN;
p = min(p,1-p);
