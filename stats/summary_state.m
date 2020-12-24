clear
clc
close all

hind = 0;

timelims = [0 5];
speclims = [.5 .92];

h = zeros(20,1);
hfig = h;

bootN = 5000;

%TIMING VAR. 

load('timeanl_DIR.mat')
d = (sqrt(paramat(:,1))-sqrt(paramat(:,2)))./sqrt(paramat(:,1));
dhatdist = bootstrp_grouped_dist(d,birdall,bootN,@median);
p = sum(dhatdist<0)/bootN;
p_psi_time_dir = min(p,1-p);
med_psi_time_dir = mean(dhatdist);
sem_psi_time_dir = std(dhatdist);
birdstats_psi_time_dir = grpmedian(d,birdall);

d = (paramat(:,5)-paramat(:,6))./paramat(:,5);
dhatdist = bootstrp_grouped_dist(d,birdall,bootN,@median);
p = sum(dhatdist<0)/bootN;
p_W_time_dir = min(p,1-p);
med_W_time_dir = mean(dhatdist);
sem_W_time_dir = std(dhatdist);
birdstats_W_time_dir = grpmedian(d,birdall);

d = (paramat(:,3)-paramat(:,4))./paramat(:,3);
dhatdist = bootstrp_grouped_dist(d,birdall,bootN,@median);
p = sum(dhatdist<0)/bootN;
p_mn_time_dir = min(p,1-p);
med_mn_time_dir = mean(dhatdist);
sem_mn_time_dir = std(dhatdist);
birdstats_mn_time_dir = grpmedian(d,birdall);

figure
hind = hind+1;

plotbird(sqrt(paramat(:,1)),sqrt(paramat(:,2)),birdall,gca)
%gscatter(sqrt(paramat(:,1)),sqrt(paramat(:,2)),birdall,'k','o',[],'off');
h(hind) = gca;
hfig(hind) = gcf;
hold on
plot(timelims,timelims,'k--','linewidth',2)
xlim(timelims)
ylim(timelims)


figure
%gscatter(paramat(:,5),paramat(:,6),birdall,'k','o',[],'off');

plotbird(paramat(:,5),paramat(:,6),birdall,gca)

hold on
plot(timelims,timelims,'k--','linewidth',2)
xlim(timelims)
ylim(timelims)

hind = hind+1;
h(hind) = gca;
hfig(hind) = gcf;

clear paramat

load('timeanl_NE.mat')
d = (sqrt(paramat(:,1))-sqrt(paramat(:,2)))./sqrt(paramat(:,1));
dhatdist = bootstrp_grouped_dist(d,birdall,bootN,@median);
p = sum(dhatdist<0)/bootN;
p_psi_time_NE = min(p,1-p);
med_psi_time_NE = mean(dhatdist);
sem_psi_time_NE = std(dhatdist);
birdstats_psi_time_NE = grpmedian(d,birdall);

d = (paramat(:,5)-paramat(:,6))./paramat(:,5);
dhatdist = bootstrp_grouped_dist(d,birdall,bootN,@median);
p = sum(dhatdist<0)/bootN;
p_W_time_NE = min(p,1-p);
med_W_time_NE = mean(dhatdist);
sem_W_time_NE = std(dhatdist);
birdstats_W_time_NE = grpmedian(d,birdall);

d = (paramat(:,3)-paramat(:,4))./paramat(:,3);
dhatdist = bootstrp_grouped_dist(d,birdall,bootN,@median);
p = sum(dhatdist<0)/bootN;
p_mn_time_NE = min(p,1-p);
med_mn_time_NE = mean(dhatdist);
sem_mn_time_NE = std(dhatdist);
birdstats_mn_time_NE = grpmedian(d,birdall);


figure
%gscatter(sqrt(paramat(:,1)),sqrt(paramat(:,2)),birdall,'k','o',[],'off');

plotbird(sqrt(paramat(:,1)),sqrt(paramat(:,2)),birdall,gca)
hold on
plot(timelims,timelims,'k--','linewidth',2)
xlim(timelims)
ylim(timelims)

hind = hind+1;
h(hind) = gca;
hfig(hind) = gcf;


figure
%gscatter(paramat(:,5),paramat(:,6),birdall,'k','o',[],'off');

plotbird(paramat(:,5),paramat(:,6),birdall,gca)

hold on
plot(timelims,timelims,'k--','linewidth',2)
xlim(timelims)
ylim(timelims)
hind = hind+1;
h(hind) = gca;
hfig(hind) = gcf;

load('timeanl_PHE.mat')
d = (sqrt(paramat(:,1))-sqrt(paramat(:,2)))./sqrt(paramat(:,1));
dhatdist = bootstrp_grouped_dist(d,birdall,bootN,@median);
p = sum(dhatdist<0)/bootN;
p_psi_time_PHE = min(p,1-p);
med_psi_time_PHE = mean(dhatdist);
sem_psi_time_PHE = std(dhatdist);
birdstats_psi_time_PHE = grpmedian(d,birdall);

d = (paramat(:,5)-paramat(:,6))./paramat(:,5);
dhatdist = bootstrp_grouped_dist(d,birdall,bootN,@median);
p = sum(dhatdist<0)/bootN;
p_W_time_PHE = min(p,1-p);
med_W_time_PHE = mean(dhatdist);
sem_W_time_PHE = std(dhatdist);
birdstats_W_time_PHE = grpmedian(d,birdall);

d = (paramat(:,3)-paramat(:,4))./paramat(:,3);
dhatdist = bootstrp_grouped_dist(d,birdall,bootN,@median);
p = sum(dhatdist<0)/bootN;
p_mn_time_PHE = min(p,1-p);
med_mn_time_PHE = mean(dhatdist);
sem_mn_time_PHE = std(dhatdist);
birdstats_mn_time_PHE = grpmedian(d,birdall);



figure
%gscatter(sqrt(paramat(:,1)),sqrt(paramat(:,2)),birdall,'k','o',[],'off');

plotbird(sqrt(paramat(:,1)),sqrt(paramat(:,2)),birdall,gca)

hold on
plot(timelims,timelims,'k--','linewidth',2)
xlim(timelims)
ylim(timelims)
hind = hind+1;
h(hind) = gca;
hfig(hind) = gcf;

figure
%gscatter(paramat(:,5),paramat(:,6),birdall,'k','o',[],'off');

plotbird(paramat(:,5),paramat(:,6),birdall,gca)

hold on
plot(timelims,timelims,'k--','linewidth',2)
xlim(timelims)
ylim(timelims)
hind = hind+1;
h(hind) = gca;
hfig(hind) = gcf;

load('timeanl_LCstim.mat')
d = (sqrt(paramat(:,1))-sqrt(paramat(:,2)))./sqrt(paramat(:,1));
dhatdist = bootstrp_grouped_dist(d,birdall,bootN,@median);
p = sum(dhatdist<0)/bootN;
p_psi_time_LC = min(p,1-p);
med_psi_time_LC = mean(dhatdist);
sem_psi_time_LC = std(dhatdist);
birdstats_psi_time_LC = grpmedian(d,birdall);

d = (paramat(:,5)-paramat(:,6))./paramat(:,5);
dhatdist = bootstrp_grouped_dist(d,birdall,bootN,@median);
p = sum(dhatdist<0)/bootN;
p_W_time_LC = min(p,1-p);
med_W_time_LC = mean(dhatdist);
sem_W_time_LC = std(dhatdist);
birdstats_W_time_LC = grpmedian(d,birdall);


d = (paramat(:,3)-paramat(:,4))./paramat(:,3);
dhatdist = bootstrp_grouped_dist(d,birdall,bootN,@median);
p = sum(dhatdist<0)/bootN;
p_mn_time_LC = min(p,1-p);
med_mn_time_LC = mean(dhatdist);
sem_mn_time_LC = std(dhatdist);
birdstats_mn_time_LC = grpmedian(d,birdall);


figure
% gscatter(sqrt(paramat(:,1)),sqrt(paramat(:,2)),birdall,'k','o',[],'off');

plotbird(sqrt(paramat(:,1)),sqrt(paramat(:,2)),birdall,gca)

hold on
plot(timelims,timelims,'k--','linewidth',2)
xlim(timelims)
ylim(timelims)
hind = hind+1;
h(hind) = gca;
hfig(hind) = gcf;

figure
%gscatter(paramat(:,5),paramat(:,6),birdall,'k','o',[],'off');

plotbird(paramat(:,5),paramat(:,6),birdall,gca)

hold on
plot(timelims,timelims,'k--','linewidth',2)
xlim(timelims)
ylim(timelims)
hind = hind+1;
h(hind) = gca;
hfig(hind) = gcf;

%SPEC VAR.

load('specanl_DIR.mat')
d = (paramat(:,1)-paramat(:,2))./paramat(:,1);
dhatdist = bootstrp_grouped_dist(d,birdall,bootN,@median);
p = sum(dhatdist<0)/bootN;
p_spec_dir = min(p,1-p);
med_spec_dir = mean(dhatdist);
sem_spec_dir = std(dhatdist);
birdstats_spec_dir = grpmedian(d,birdall);


figure
%gscatter(paramat(:,1),paramat(:,2),birdall,'k','o',[],'off');

plotbird(paramat(:,1),paramat(:,2),birdall,gca)


hold on
plot(speclims,speclims,'k--','linewidth',2)
xlim(speclims)
ylim(speclims)
hind = hind+1;
h(hind) = gca;
hfig(hind) = gcf;


load('specanl_NE.mat')
d = (paramat(:,1)-paramat(:,2))./paramat(:,1);
dhatdist = bootstrp_grouped_dist(d,birdall,bootN,@median);
p = sum(dhatdist<0)/bootN;
p_spec_NE = min(p,1-p);
med_spec_NE = mean(dhatdist);
sem_spec_NE = std(dhatdist);
birdstats_spec_NE = grpmedian(d,birdall);

figure
%gscatter(paramat(:,1),paramat(:,2),birdall,'k','o',[],'off');
plotbird(paramat(:,1),paramat(:,2),birdall,gca)


hold on
plot(speclims,speclims,'k--','linewidth',2)
xlim(speclims)
ylim(speclims)
hind = hind+1;
h(hind) = gca;
hfig(hind) = gcf;



load('specanl_PHE.mat')
d = (paramat(:,1)-paramat(:,2))./paramat(:,1);
dhatdist = bootstrp_grouped_dist(d,birdall,bootN,@median);
p = sum(dhatdist<0)/bootN;
p_spec_PHE = min(p,1-p);
med_spec_PHE = mean(dhatdist);
sem_spec_PHE = std(dhatdist);
birdstats_spec_PHE = grpmedian(d,birdall);

figure
%gscatter(paramat(:,1),paramat(:,2),birdall,'k','o',[],'off');
plotbird(paramat(:,1),paramat(:,2),birdall,gca)

hold on
plot(speclims,speclims,'k--','linewidth',2)
xlim(speclims)
ylim(speclims)
hind = hind+1;
h(hind) = gca;
hfig(hind) = gcf;


load('specanl_LCstim.mat')
d = (paramat(:,1)-paramat(:,2))./paramat(:,1);
dhatdist = bootstrp_grouped_dist(d,birdall,bootN,@median);
p = sum(dhatdist<0)/bootN;
p_spec_LC = min(p,1-p);
med_spec_LC = mean(dhatdist);
sem_spec_LC = std(dhatdist);
birdstats_spec_LC = grpmedian(d,birdall);

figure
% gscatter(paramat(:,1),paramat(:,2),birdall,'k','o',[],'off');

plotbird(paramat(:,1),paramat(:,2),birdall,gca);

hold on
plot(speclims,speclims,'k--','linewidth',2)
xlim(speclims)
ylim(speclims)
hind = hind+1;
h(hind) = gca;
hfig(hind) = gcf;


h = h(find(h));
hfig = hfig(find(hfig));

for hind = 1:numel(hfig)
    set(hfig(hind),'Position',[ 440   691   150   107])
    set(h(hind),'Fontsize',14,'box','off','ticklength',[0 0])
end