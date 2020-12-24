clear
clc
close all

% conditions: 
% 1: undir. no stim  
% 2: undir. stim.
% 3: undir. sal
% 4: dir. sal
% 5: undir. NE
% 6: dir. PHE

% measure type:
% 1: spec var.
% 2: spec mn.
% 3: dur mn.
% 4: dur ind var.
% 5: dur W var.

birdmain = [];
intmain = [];
measuremain = [];
measuretype = [];
condtype = [];

load('specanl_LCstim.mat')

birdmain = [birdmain;birdarr];
intmain = [intmain;sylarr];

load('specanl_dir.mat')

load('specanl_NE.mat')
load('specanl_PHE.mat')