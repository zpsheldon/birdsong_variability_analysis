clear
clc
close all

load /Users/cg/Documents/MATLAB/NE_data/OR_46/specdata_honed3.mat
load /Users/cg/Documents/MATLAB/NE_data/OR_13/OR_13_F_baseline_all/optsample.mat propstrct

freqs = propstrct.fs*[0:propstrct.f_winlen/2]/propstrct.f_winlen;
freqinds = find(freqs >= propstrct.freqmin & freqs <= propstrct.freqmax);
freqs = freqs(freqinds);

sylind = 4;
sylpts = 5:7;

sylmat_NE = specstrct_F_sal.spalgnarr{sylind};
sylmat_NE = sylmat_NE(:,:,sylpts);

sylmat_sal = specstrct_U_sal.spalgnarr{sylind};
sylmat_sal = sylmat_sal(:,:,sylpts);

sylmn_NE = specstrct_F_sal.spmnarr{sylind}(:,sylpts);
sylmn_sal = specstrct_U_sal.spmnarr{sylind}(:,sylpts);

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

indsval_sal = FF_sal>2800;
indsval_NE = FF_NE>2800;
FF_sal = FF_sal(indsval_sal);
FF_NE = FF_NE(indsval_NE);
dist_sal = dist_sal(indsval_sal);
dist_NE = dist_NE(indsval_NE);

figure; plot(abs(FF_sal-mean(FF_sal)),dist_sal,'ko')
hold on
plot(abs(FF_NE-mean(FF_NE)),dist_NE,'o','color',[.2 .7 .2])


% xlabel('deviation in fundamental frequency (Hz)','fontsize',14)
% ylabel('general spectro-temporal distance (A.U.)','fontsize',14)
set(gca,'box','off','ticklength',[0 0],'fontsize',12)

legend({'alone','female present'},'fontsize',12)



xvals = 2800:20:3000;

figure
subplot(2,1,1)
[N,xout] = hist(FF_sal,xvals);
bar(xout,N,'facecolor','k','barwidth',.8)

subplot(2,1,2)
[N,xout] = hist(FF_NE,xvals);

bar(xout,N,'facecolor','k','barwidth',.8)

for i = 1:2
    subplot(2,1,i)
    set(gca,'ticklength',[0 0],'box','off','fontsize',12)
    xlim([xvals(1) xvals(end)])
end


set(gcf,'Position',[440   642   148   156])

figure
imagesc(specstrct_U_sal.spmnarr{sylind})
set(gca,'box','off','ticklength',[0 0],'ydir','normal','xtick',[],'ytick',[])
set(gcf,'Position',[440   691   171   107])


% figure; imagesc(specstrct_U_NE.spmnarr{5})