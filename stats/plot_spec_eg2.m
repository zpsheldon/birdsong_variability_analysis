clear
clc
close all

load /Users/cg/Documents/MATLAB/NE_data/RD_2/specdata.mat
load /Users/cg/Documents/MATLAB/NE_data/RD_2/STIM/cliptraining.mat propstrct

freqs = propstrct.fs*[0:propstrct.f_winlen/2]/propstrct.f_winlen;
freqinds = find(freqs >= propstrct.freqmin & freqs <= propstrct.freqmax);
freqs = freqs(freqinds);

sylind = 6;
sylpts = 2:7;

sylmat_NE = specstrct_U_stim.spalgnarr{sylind};
sylmat_NE = sylmat_NE(:,:,sylpts);

sylmat_sal = specstrct_U_nostim.spalgnarr{sylind};
sylmat_sal = sylmat_sal(:,:,sylpts);

sylmn_NE = specstrct_U_stim.spmnarr{sylind}(:,sylpts);
sylmn_sal = specstrct_U_nostim.spmnarr{sylind}(:,sylpts);

N_NE = size(sylmat_NE,1);
N_sal = size(sylmat_sal,1);

FF_NE = zeros(N_NE,1);
FF_sal = zeros(N_sal,1);

dist_NE = FF_NE;
dist_sal = FF_sal;

indsval = 1:numel(freqs);

for n = 1:N_sal
   syltmp = squeeze(sylmat_sal(n,:,:));  
   [fund,pwrat,ceps] = fund_freq(syltmp(indsval,:),freqs(indsval));
   FF_sal(n) = mean(fund);
   dist_sal(n) = sqrt(mean((syltmp(:)-sylmn_sal(:)).^2)); 
end

for n = 1:N_NE
   syltmp = squeeze(sylmat_NE(n,:,:));  
   [fund,pwrat,ceps] = fund_freq(syltmp(indsval,:),freqs(indsval));
   FF_NE(n) = mean(fund);
   dist_NE(n) = sqrt(mean((syltmp(:)-sylmn_NE(:)).^2)); 
end

figure; plot(abs(FF_sal-mean(FF_sal)),dist_sal,'ko')
% hold on
% plot(abs(FF_NE-mean(FF_NE)),dist_NE,'o')

% xlabel('deviation in fundamental frequency (Hz)','fontsize',14)
% ylabel('general spectro-temporal distance (A.U.)','fontsize',14)
set(gca,'box','off','ticklength',[0 0],'fontsize',12)

set(gcf,'Position',[440   691   171   107])

figure
imagesc(specstrct_U_nostim.spmnarr{sylind})
set(gca,'box','off','ticklength',[0 0],'ydir','normal','xtick',[],'ytick',[])
set(gcf,'Position',[440   691   171   107])


% legend({'saline','NE'},'fontsize',12)

% figure; imagesc(specstrct_U_NE.spmnarr{5})