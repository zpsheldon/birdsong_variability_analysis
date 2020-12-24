clear
clc
close all

% BR2: have data and anl
% OR13: have data and anl
% BR26: have data, NOT anl
% OR2: have data, NOT anl
% WH57: have data NOT anl
% Y437: no data
% Y439: no data

maindir = ['/Users/cg/Documents/MATLAB/NE_data/'];

% wavdirs = {'BR26/UNE','BR26/USAL','OR2/USAL','OR2/UNE','WH57/USAL','WH57/UNE'};

wavdirs = {'BL_16/STIM','BL_16/NOSTIM','SI_026/STIM','SI_026/NOSTIM','PU_31/STIM','PU_31/NOSTIM','RD_2/STIM','RD_2/NOSTIM'};

for wavdirind = 1:numel(wavdirs)
    
    wavdir = [maindir wavdirs{wavdirind}];
    dirstrct = mk_wavdirinfo(wavdir,'Timestamps');
    clipstrct = mk_clipdirinfo(wavdir,'Timestamps');
    
end
