clear
clc
close all

birdarr = {'OR_46','BR26','WH57','BR_2','WH27'};

birdind = 5;

load(['./' birdarr{birdind} '/specdata_honed.mat'])
%load(['./' birdarr{birdind} '/specdata.mat'])

[stdtmp_exp,distmp_exp,stdtot_exp,spec_exp,specdev_exp,timeinds_exp] = spec_regression(specstrct_F_sal.spalgnarr);
% 
% indsrand = randperm(size(specstrct_U_sal.spalgnarr{1},1));
% indsrand = indsrand(1:20);

% for i = 1:numel(specstrct_U_sal.spalgnarr)
%     specstrct_U_sal.spalgnarr{i} = specstrct_U_sal.spalgnarr{i}(indsrand,:,:);
% end

[stdtmp_sal,distmp_sal,stdtot_sal,spec_sal,specdev_sal,timeinds_sal] = spec_regression(specstrct_F_PHE.spalgnarr);

[seqN_sal,freqN_sal,timeN_sal] = size(specdev_sal);
[seqN_exp,freqN_exp,timeN_exp] = size(specdev_exp);

freqs = 500+(9000-500)*[0:freqN_sal]/freqN_sal;

ffmat_exp = zeros(seqN_exp,timeN_exp);
ffmat_sal = zeros(seqN_sal,timeN_sal);
pwrmat_exp = ffmat_exp;
pwrmat_sal = ffmat_sal;

for seqind = 1:seqN_exp
    spec = squeeze(spec_exp(seqind,:,:));
    featmat = specfeat2(spec,freqs);
    ffmat_exp(seqind,:) = featmat(1,:);
    pwrmat_exp(seqind,:) = featmat(3,:);
end

for seqind = 1:seqN_sal
    spec = squeeze(spec_sal(seqind,:,:));
    featmat = specfeat2(spec,freqs);
    ffmat_sal(seqind,:) = featmat(1,:);
    pwrmat_sal(seqind,:) = featmat(3,:);
end

% vartmp_sal = mean(squeeze(var(specdev_sal,1,1)));
% vartmp_exp = mean(squeeze(var(specdev_exp,1,1)));

vartmp_sal = mean(squeeze(var(specdev_sal,1,1)));
vartmp_exp = mean(squeeze(var(specdev_exp,1,1)));

spmn_sal = squeeze(mean(spec_sal,1));
spmn_exp = squeeze(mean(spec_exp,1));

boundaries = zeros(1,numel(timeinds_sal)+1);
boundaries(1) = 1;

for sylind = 1:numel(timeinds_sal)
    boundaries(sylind+1) = timeinds_sal{sylind}(end); 
end

t = [1:size(spmn_sal,2)];

figure
subplot(5,1,1)
imagesc(t,freqs,spmn_sal,[0 5])
set(gca,'ydir','normal')

title('PHE mean')

subplot(5,1,2)
imagesc(t,freqs,spmn_exp,[0 5])
set(gca,'ydir','normal')

title('sal mean')

subplot(5,1,3)
plot(vartmp_sal)
hold on
plot(vartmp_exp)
xlim([1 size(spmn_exp,2)])

legend({'PHE','sal'})

subplot(5,1,4)
plot(mean(ffmat_sal))
hold on
plot(mean(ffmat_exp))
xlim([1 size(spmn_exp,2)])

subplot(5,1,5)
plot(std(ffmat_sal))
hold on
plot(std(ffmat_exp))
xlim([1 size(spmn_exp,2)])


for i = 1:5
   subplot(5,1,i)
   ylims = get(gca,'ylim');
   hold on
   plot([boundaries;boundaries],repmat(ylims(:),1,numel(boundaries)),'k-','linewidth',2)
    
end

set(gcf,'Position',[180   385   560   420])

bintarg = 22;

bintarg = 10;

timeslice_sal = squeeze(spec_sal(:,:,bintarg));
timeslice_exp = squeeze(spec_exp(:,:,bintarg));
% 
% figure
% subplot(4,1,1)
% imagesc(timeslice_sal,[0 5])
% title('undir')
% 
% subplot(4,1,2)
% imagesc(timeslice_exp,[0 5])
% title('dir')
% 
% subplot(4,1,3)
% plot(var(timeslice_sal))
% hold on
% plot(var(timeslice_exp))
% legend({'undir','dir'})
% title('variance')
% 
% subplot(4,1,4)
% plot(mean(timeslice_sal))
% hold on
% plot(mean(timeslice_exp))
% legend({'undir','dir'})
% 
% title('mean')
% 
% set(gcf,'Position',[752   385   560   420])
% 
% %specdevsal_coll = reshape(specdev_sal,seqN_sal,freqN_sal*timeN_sal);