function config_strct = ABconfig

config_strct.freqwin = 10;
config_strct.freqmax = 9000;
config_strct.freqmin = 500;
config_strct.timewin = 5;

config_strct.amp_cutoff = 0.025;

% config_strct.amp_cutoff = 0.0025;

config_strct.ampwin = 1.5;

config_strct.gap_thresh = 15;
config_strct.amp_thresh = .5;

%config_strct.amp_thresh = .2;

config_strct.logopt = 1;
config_strct.normopt = 1;

config_strct.clip_amp_cutoff = 0.01;
config_strct.clip_minlen = 10;
config_strct.clip_maxlen = 500;

config_strct.filt_order = 5;

config_strct.freqwin_fn = 5;
config_strct.timewin_fn = 0.16;
config_strct.smthwn = 10;


% ORIGINAL PARAMETERS FROM SCHMIDT LAB:
% 
% config_strct.freqwin = 10;
% config_strct.freqmax = 9000;
% config_strct.freqmin = 500;
% config_strct.timewin = 5;
% 
% config_strct.amp_cutoff = 0.025;
% 
% config_strct.ampwin = 1.5;
% 
% config_strct.gap_thresh = 7;
% config_strct.amp_thresh = .5;
% 
% config_strct.logopt = 1;
% config_strct.normopt = 1;
% 
% config_strct.clip_amp_cutoff = 0.01;
% config_strct.clip_minlen = 10;
% config_strct.clip_maxlen = 500;
% 
% config_strct.filt_order = 5;
% 
% config_strct.freqwin_fn = 5;
% config_strct.timewin_fn = 0.16;
% config_strct.smthwn = 10;