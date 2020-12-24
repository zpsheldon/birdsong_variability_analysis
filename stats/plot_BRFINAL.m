clear
clc
close all

load('timeanl_NE.mat')

birdinds = find(strcmp(birdall,'BRFINAL'));
params = paramat(birdinds,:);
syls = sylall(birdinds);

figure
plot(params(:,1),'-o')
hold on
plot(params(:,2),'-+')