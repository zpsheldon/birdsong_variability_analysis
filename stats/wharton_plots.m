clear
clc
close all

load songFreqData_LCstim.mat



figure; 

bar([1 2],[median(sngRateUndir) median(sngRateDir)],'barwidth',.8,'facecolor','k')
set(gca,'box','off','ticklength',[0 0],'fontsize',14,'xtick',[1 2],'xticklabel',{'no NE stim','NE stim'})
ylabel('songs/minute','fontsize',16)

set(gcf,'Position',[440   500   271   200])


load('specanl_LCstim.mat')

figure; 

bar([1 2],[median(paramat(:,1:2))],'barwidth',.8,'facecolor','k')
set(gca,'box','off','ticklength',[0 0],'fontsize',14,'xtick',[1 2],'xticklabel',{'no NE stim','NE stim'})
ylabel('spectral var.','fontsize',16)
ylim([.7 .78])

set(gcf,'Position',[440   500   271   200])


load('timeanl_LCstim.mat')

figure; 

bar([1 2],[median(paramat(:,1:2))],'barwidth',.8,'facecolor','k')
set(gca,'box','off','ticklength',[0 0],'fontsize',14,'xtick',[1 2],'xticklabel',{'no NE stim','NE stim'})
ylabel('rhythmic var.','fontsize',16)
ylim([.2 1.2])

set(gcf,'Position',[440   500   271   200])

