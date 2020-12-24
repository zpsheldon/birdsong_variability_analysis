#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 28 13:56:57 2019

@author: Zach Sheldon
"""

# plotting functions
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from scipy.stats import ttest_rel
from scipy.stats import sem

################################## RAW DATA #############################################

# song frequency - BR2, BR26, OR2, OR13, OR46, WH27, WH57

# LC Stimulation vs. No Stimulation

# Undirected - - RD02 (old), PU31 (old), BL16 (old), SI26 (old), RD08 (new), WH09 (new)
lc_nostim_undir_rd02 = [0.6494, 0.620477886643942, 0.548859944421228, 0.671251315744692, 0.628774789380207, 0.659000758745500]
lc_stim_undir_rd02 = [0.639369741800785, 0.595635411764726, 0.529020762363586, 0.654010831140285, 0.620755754180120, 0.650853166390532]

lc_nostim_undir_pu31 = [0.615563251860123, 0.678034748840814	, 0.659994768702187, 0.618504643294967, 0.650795253230500, 0.644207049646824, 0.662579264122903]
lc_stim_undir_pu31 = [0.599167128262597, 0.655138263570173, 0.653226995883096, 0.619839986372945, 0.639964506682945, 0.636540309176339, 0.648634370701892]

lc_nostim_undir_bl16 = [0.640639338380572, 0.616520723238547, 0.593281552234197, 0.661168895973408, 0.818140543281690]
lc_stim_undir_bl16 = [0.647936945039651, 0.616606897921562, 0.594290391360265, 0.643902927259642, 0.774782930628622]

lc_nostim_undir_si26 = [0.589984412458911, 0.665677524840571, 0.679033536637287, 0.612314111869927, 0.740809561237048, 0.707965713178863, 0.849228457676734]
lc_stim_undir_si26 = [0.56097114096733, 0.651099531386898, 0.665060267966882, 0.572585225199502, 0.691410767572151, 0.615149274512045, 0.782965970441265]

lc_nostim_undir_rd08 = [0.6719, 0.5834, 0.5311, 0.7262]
lc_stim_undir_rd08 = [0.6320, 0.5507, 0.5207, 0.6729]

lc_nostim_undir_wh09 = [0.6615, 0.5309, 0.5854, 0.5516]
lc_stim_undir_wh09 = [0.6960, 0.5697, 0.7112, 0.6411]

lc_stim_undir_overall = lc_stim_undir_rd02 + lc_stim_undir_pu31 + lc_stim_undir_bl16 + lc_stim_undir_si26
lc_nostim_undir_overall = lc_nostim_undir_rd02 + lc_nostim_undir_pu31 + lc_nostim_undir_bl16 + lc_nostim_undir_si26
lc_stim_undir_miss_overall = lc_stim_undir_wh09
lc_nostim_undir_miss_overall = lc_nostim_undir_wh09

# Directed - RD08 (LC hit), WH09 (LC miss)
lc_nostim_dir_rd08 = [0.6311, 0.5632, 0.4776, 0.6547]
lc_stim_dir_rd08 = [0.6882, 0.6021, 0.5643, 0.6976]

lc_nostim_dir_wh09 = [0.6905, 0.5853, 0.6283, 0.5442]
lc_stim_dir_wh09 = [0.6348, 0.5454, 0.6085, 0.4517]

lc_stim_dir_overall = lc_stim_dir_rd08
lc_nostim_dir_overall = lc_nostim_dir_rd08
lc_stim_dir_miss_overall = lc_stim_dir_wh09
lc_nostim_dir_miss_overall = lc_nostim_dir_wh09

# NE Infusion vs. Saline

# Undirected - BR_2, OR_13, BR26, OR2, WH57, Y437, WH27, BR0
ne_sal_undir_br02 = [0.647946405268450, 0.616162601182828, 0.663480524399165, 0.582218884150931, 0.662177584767032, 0.693383742418067]
ne_inf_undir_br02 = [0.577299373121745, 0.559785740719010, 0.595498259407118, 0.539754123042061, 0.632090394084641, 0.603001055813073]

ne_inf_undir_br0 = [0.616773829502707, 0.615749017198440, 0.602455833504472, 0.629750015166038, 0.652666427593210]
ne_sal_undir_br0 = [0.630242078140791, 0.628769922754563, 0.651758630349732, 0.669034945889631, 0.706120134779261]

ne_inf_undir_or13 = [0.599388749031826, 0.665721644290975, 0.594411813706152, 0.580233248715771]
ne_sal_undir_or13 = [0.617134242186315, 0.689862787559255, 0.663679960430387, 0.658618955156036]

ne_inf_undir_br26 = [0.679755544050973, 0.669075855736476]
ne_sal_undir_br26 = [0.665538759394075, 0.676824745831506]

ne_inf_undir_or02 = [0.650113105722325, 0.624695087179132, 0.648535887910220, 0.687622304607323]
ne_sal_undir_or02 = [0.674810163639679, 0.642769663140988, 0.680519869465259, 0.715878608540469]

ne_inf_undir_wh27 = [0.633576317889146, 0.633718012975067, 0.671235872438273]
ne_sal_undir_wh27 = [0.716123097537054, 0.683453292268733, 0.773143372436969]

ne_inf_undir_wh57 = [0.586562302805536, 0.700287927605865, 0.629168844646946, 0.686009280113732, 0.744108339343115]
ne_sal_undir_wh57 = [0.651583588184139, 0.689030426588848, 0.643299218932953, 0.708329181737318, 0.713112953768411]

ne_inf_undir_y437 = [0.711130877354350, 0.746263430687309, 0.776973398260017, 0.751116651671179, 0.748014866644673, 0.746081875384835]
ne_sal_undir_y437 = [0.714571353574320, 0.735238327755120, 0.776731647608070, 0.761119148081544, 0.753573968377751, 0.748643145527786]

ne_inf_undir_overall = ne_inf_undir_br02 + ne_inf_undir_br0 + ne_inf_undir_or13 + ne_inf_undir_br26 + ne_inf_undir_or02 + ne_inf_undir_wh27 + ne_inf_undir_wh57 + ne_inf_undir_y437
ne_sal_undir_overall = ne_sal_undir_br02 + ne_sal_undir_br0 + ne_sal_undir_or13 + ne_sal_undir_br26 + ne_sal_undir_or02 + ne_sal_undir_wh27 + ne_sal_undir_wh57 + ne_sal_undir_y437

# Directed - BR_2, OR_13, BR26, OR2, WH27, WH57
ne_sal_dir_br02 = [0.614642786615478, 0.606367769383930, 0.652151545139075, 0.569817999157838, 0.674621401691125, 0.702524925602959]
ne_inf_dir_br02 = [0.590319003170878, 0.547701332555572, 0.563399618495704, 0.515593855312563, 0.552361305792287, 0.531508934303613]

ne_sal_dir_br26 = [0.667222929798158, 0.681040387361729]
ne_inf_dir_br26 = [0.652426584741179, 0.661459136685726]

ne_sal_dir_or13 = [0.612605742404265, 0.704704652871119, 0.637067368049386, 0.637947912065505]
ne_inf_dir_or13 = [0.600270477795239, 0.613470850440367, 0.576974573182433, 0.567601959248593]

ne_sal_dir_or02 = [0.660801413064616, 0.638573811488045, 0.645981770930920, 0.673077238783320]
ne_inf_dir_or02 = [0.633599138739974, 0.615257301254814, 0.613916885180071, 0.689270632702252]

ne_sal_dir_wh27 = [0.703113965756233, 0.670391541710256, 0.753243858463428]
ne_inf_dir_wh27 = [0.730530974574548, 0.718533616804567, 0.773176793244452]

ne_sal_dir_wh57 = [0.650956691341111, 0.688692907062697, 0.643224210805544, 0.708634407242040, 0.712659715795996]
ne_inf_dir_wh57 = [0.658962243300032, 0.684297372869443, 0.632364651820652, 0.718817488170459, 0.711379762119059]

ne_inf_dir_overall = ne_inf_dir_br26 + ne_inf_dir_or13 + ne_inf_dir_or02 + ne_inf_dir_wh27 + ne_inf_dir_wh57
ne_sal_dir_overall = ne_sal_dir_br26 + ne_sal_dir_or13 + ne_sal_dir_or02 + ne_sal_dir_wh27 + ne_sal_dir_wh57

# Phentolamine Infusion vs. Saline


# Directed - OR46, BR26, WH27, WH57
phe_sal_dir_or46 = [0.523252268551732, 0.491019748614369, 0.559294571029945, 0.599779093339930, 0.530967992642175, 0.435030361928949, 0.527434441680192, 0.457547858753285]
phe_inf_dir_or46 = [0.574288289612474, 0.530726626942409, 0.609735776879352, 0.625028338324387, 0.575323401088860, 0.474454549581341, 0.630429090315022, 0.527531027290634]

phe_sal_dir_br26 = [0.652426584741179, 0.661459136685726]
phe_inf_dir_br26 = [0.662537002055870, 0.667327975899367]

phe_sal_dir_wh27 = [0.730530974574548, 0.718533616804567, 0.773176793244452]
phe_inf_dir_wh27 = [0.723509620911941, 0.717163180337522, 0.753717562470245]

phe_sal_dir_wh57 = [0.658962243300032, 0.684297372869443, 0.632364651820652, 0.718817488170459, 0.711379762119059]
phe_inf_dir_wh57 = [0.685652236241114, 0.770066374778060, 0.805444389439339, 0.817972664177677, 0.846778178726015]

phe_inf_overall = phe_inf_dir_or46 + phe_inf_dir_br26 + phe_inf_dir_wh27 + phe_inf_dir_wh57
phe_sal_overall = phe_sal_dir_or46 + phe_sal_dir_br26 + phe_sal_dir_wh27 + phe_sal_dir_wh57

# Song Frequency

# LC - RD02, PU31, SI026, BL16
rd08_numsongs_undir_lcstim = (104+69+324)/3.0
rd08_numsongs_undir_nostim = (126+26+26)/3.0
rd08_songfreq_undir_lcstim = rd08_numsongs_undir_lcstim/3.0/60.0
rd08_songfreq_undir_nostim = rd08_numsongs_undir_nostim/3.0/60.0

wh09_numsongs_undir_lcstim = 216/3.0
wh09_numsongs_undir_nostim = 350/3.0
wh09_songfreq_undir_lcstim = wh09_numsongs_undir_lcstim/3.0/60.0
wh09_songfreq_undir_nostim = wh09_numsongs_undir_nostim/3.0/60.0

song_freq_undir_lcstim = [0.918666124585571, 2.02833448478658, 0.865384615644225, 2.14881678417931, rd08_songfreq_undir_lcstim]
song_freq_undir_nostim = [0.321274650679320, 0.686962486149005, 0.568163516250445, 0.402662130735009, rd08_songfreq_undir_nostim]
rd02_songfreq_undir_lcstim = 0.918666124585571
rd02_songfreq_undir_nostim = 0.321274650679320
pu31_songfreq_undir_lcstim = 2.02833448478658
pu31_songfreq_undir_nostim = 0.686962486149005
si026_songfreq_undir_lcstim = 0.865384615644225
si026_songfreq_undir_nostim = 0.568163516250445
bl16_songfreq_undir_lcstim = 2.14881678417931
bl16_songfreq_undir_nostim = 0.402662130735009

songfreq_lcstim_undir_overall = [rd08_songfreq_undir_lcstim, rd02_songfreq_undir_lcstim, bl16_songfreq_undir_lcstim, pu31_songfreq_undir_lcstim, si026_songfreq_undir_lcstim]
songfreq_lcstim_undir_sem = np.std(songfreq_lcstim_undir_overall) / np.sqrt(5)

songfreq_nostim_undir_overall = [rd08_songfreq_undir_nostim, rd02_songfreq_undir_nostim, bl16_songfreq_undir_nostim, pu31_songfreq_undir_nostim, si026_songfreq_undir_nostim]
songfreq_nostim_undir_sem = np.std(songfreq_nostim_undir_overall) / np.sqrt(5)
songfreq_lc_undir_diffs = [stim-nostim for stim,nostim in zip(songfreq_lcstim_undir_overall, songfreq_nostim_undir_overall)]
songfreq_lc_undir_mean_diff = np.mean(songfreq_lc_undir_diffs)
songfreq_lc_undir_std_diff = np.std(songfreq_lc_undir_diffs)
songfreq_lc_undir_sem_diff = np.std(songfreq_lc_undir_diffs) / np.sqrt(len(songfreq_lc_undir_diffs))

songfreq_lc_undir_t_val, songfreq_lc_undir_p_val = ttest_rel(songfreq_lcstim_undir_overall,songfreq_nostim_undir_overall)
songfreq_lc_undir_p_val = songfreq_lc_undir_p_val / 2.0

# NE -
song_freq_undir_ne = []
song_freq_undir_sal = []
br02_songfreq_undir_ne = 0.807983648457375
br02_songfreq_undir_sal = 0.408428344457098
br0_songfreq_undir_ne = 0.969503870033304
br0_songfreq_undir_sal = 0.212379742251317
br26_songfreq_undir_ne = 1.65271862657103
br26_songfreq_undir_sal = 2.06256640367153
or02_songfreq_undir_ne = 0.386895938339382
or02_songfreq_undir_sal = 0.422100205880415
or13_songfreq_undir_ne = 0.510062456885792
or13_songfreq_undir_sal = 1.14080793496730
wh27_songfreq_undir_ne = 1.53783063348285
wh27_songfreq_undir_sal = 1.86640471960032
wh57_songfreq_undir_ne = 0.155933171505229
wh57_songfreq_undir_sal = 0.273541114124578
y437_songfreq_undir_ne = 1.68961033139108
y437_songfreq_undir_sal = 1.14292275909326

songfreq_ne_undir_overall = [br02_songfreq_undir_ne, br0_songfreq_undir_ne, br26_songfreq_undir_ne, or02_songfreq_undir_ne, or13_songfreq_undir_ne, wh27_songfreq_undir_ne, wh57_songfreq_undir_ne, y437_songfreq_undir_ne]
songfreq_sal_undir_overall = [br02_songfreq_undir_sal, br0_songfreq_undir_sal, br26_songfreq_undir_sal, or02_songfreq_undir_sal, or13_songfreq_undir_sal, wh27_songfreq_undir_sal, wh57_songfreq_undir_sal, y437_songfreq_undir_sal]

songfreq_ne_undir_mean = np.mean([br02_songfreq_undir_ne, br0_songfreq_undir_ne, br26_songfreq_undir_ne, or02_songfreq_undir_ne, or13_songfreq_undir_ne, wh27_songfreq_undir_ne, wh57_songfreq_undir_ne, y437_songfreq_undir_ne])
songfreq_sal_undir_mean = np.mean([br02_songfreq_undir_sal, br0_songfreq_undir_sal, br26_songfreq_undir_sal, or02_songfreq_undir_sal, or13_songfreq_undir_sal, wh27_songfreq_undir_sal, wh57_songfreq_undir_sal, y437_songfreq_undir_sal])
songfreq_ne_undir_sem = np.std([br02_songfreq_undir_ne, br0_songfreq_undir_ne, br26_songfreq_undir_ne, or02_songfreq_undir_ne, or13_songfreq_undir_ne, wh27_songfreq_undir_ne, wh57_songfreq_undir_ne, y437_songfreq_undir_ne]) / np.sqrt(8)
songfreq_sal_undir_sem = np.std([br02_songfreq_undir_sal, br0_songfreq_undir_sal, br26_songfreq_undir_sal, or02_songfreq_undir_sal, or13_songfreq_undir_sal, wh27_songfreq_undir_sal, wh57_songfreq_undir_sal, y437_songfreq_undir_sal]) / np.sqrt(8)
songfreq_ne_undir_diffs = [ne-sal for ne,sal in zip(songfreq_ne_undir_overall, songfreq_sal_undir_overall)]
songfreq_ne_undir_mean_diff = np.mean(songfreq_ne_undir_diffs)
songfreq_ne_undir_std_diff = np.std(songfreq_ne_undir_diffs)
songfreq_ne_undir_sem_diff = np.std(songfreq_ne_undir_diffs) / np.sqrt(len(songfreq_ne_undir_diffs))
songfreq_ne_undir_t_val, songfreq_ne_undir_p_val = ttest_rel(songfreq_ne_undir_overall,songfreq_sal_undir_overall)
songfreq_ne_undir_p_val = songfreq_ne_undir_p_val / 2.0

# PHE -
song_freq_dir_phe = []
song_freq_dir_sal = []
or46_songfreq_dir_phe = 2.70749398327487
or46_songfreq_dir_sal = 0.920502087516287
br26_songfreq_dir_phe = 4.10646387704734
br26_songfreq_dir_sal = 4.13108838917451
wh27_songfreq_dir_phe = 5.94827588995881
wh27_songfreq_dir_sal = 1.62219101587864
wh57_songfreq_dir_phe = 0.388349516220508
wh57_songfreq_dir_sal = 1.35458167443262
songfreq_dir_phe_overall = [or46_songfreq_dir_phe, br26_songfreq_dir_phe, wh27_songfreq_dir_phe, wh57_songfreq_dir_phe]
songfreq_dir_sal_overall = [or46_songfreq_dir_sal, br26_songfreq_dir_sal, wh27_songfreq_dir_sal, wh57_songfreq_dir_sal]

songfreq_phe_dir_mean = np.mean(songfreq_dir_phe_overall)
songfreq_sal_dir_mean = np.mean(songfreq_dir_sal_overall)
songfreq_phe_sem = np.std([or46_songfreq_dir_phe, br26_songfreq_dir_phe, wh27_songfreq_dir_phe, wh57_songfreq_dir_phe])/np.sqrt(4)
songfreq_sal_sem = np.std([or46_songfreq_dir_sal, br26_songfreq_dir_sal, wh27_songfreq_dir_sal, wh57_songfreq_dir_sal])/np.sqrt(4)
songfreq_phe_dir_diffs = [phe-sal for phe,sal in zip(songfreq_dir_phe_overall, songfreq_dir_sal_overall)]
songfreq_phe_dir_mean_diff = np.mean(songfreq_phe_dir_diffs)
songfreq_phe_dir_std_diff = np.std(songfreq_phe_dir_diffs)
songfreq_phe_dir_sem_diff = np.std(songfreq_phe_dir_diffs) / np.sqrt(len(songfreq_phe_dir_diffs))

songfreq_phe_dir_t_val, songfreq_phe_dir_p_val = ttest_rel(songfreq_dir_phe_overall,songfreq_dir_sal_overall)
songfreq_phe_dir_p_val = songfreq_phe_dir_p_val / 2.0


############################# FIGURE PLOTS ###########################

# Stats

# Undirected LC Stim vs. No Stim
plt.figure(figsize=(12,8))
plt.plot(lc_nostim_undir_rd02, lc_stim_undir_rd02, 'o', label='RD02', markersize=14)
plt.plot(lc_nostim_undir_pu31, lc_stim_undir_pu31, 'o', label='PU31', markersize=14)
plt.plot(lc_nostim_undir_bl16, lc_stim_undir_bl16, 'o', label='BL16', markersize=14)
plt.plot(lc_nostim_undir_si26, lc_stim_undir_si26, 'o', label='SI026', markersize=14)
plt.plot(lc_nostim_undir_rd08, lc_stim_undir_rd08, 'o', label='RD08', markersize=14)
plt.plot(lc_nostim_undir_wh09, lc_stim_undir_wh09, 'ko', fillstyle='none', label='WH09 (LC Miss)', markersize=14)
plt.plot([0.4, 0.9], [0.4, 0.9], 'k--')
plt.xlim([0.5, 0.9])
plt.ylim([0.5, 0.9])
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.legend(fontsize=14)
plt.title('Spectral Variability (A.U.) - Undirected Song', fontsize=20)
plt.xlabel('No Stim',fontsize=20)
plt.ylabel('LC Stim',fontsize=20)
plt.tight_layout()
plt.savefig('sv_undirected_lc.png', dpi=300)

# Directed LC Stim vs. No Stim
plt.figure(figsize=(12,8))
plt.plot(lc_nostim_dir_rd08, lc_stim_dir_rd08, 'o', label='RD08 (LC Hit)', markersize=14)
plt.plot(lc_nostim_dir_wh09, lc_stim_dir_wh09, 'ko', fillstyle='none', label='WH09 (LC Miss)', markersize=14)
plt.plot([0.4, 0.9], [0.4, 0.9], 'k--')
plt.xlim([0.5, 0.9])
plt.ylim([0.5, 0.9])
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.legend(fontsize=14)
plt.title('Spectral Variability (A.U.) - Directed Song', fontsize=20)
plt.xlabel('No Stim',fontsize=20)
plt.ylabel('LC Stim',fontsize=20)
plt.savefig('sv_directed_lc.png', dpi=300)

# Undirected NE Infusion vs. Saline
plt.figure(figsize=(12,8))
plt.plot(ne_sal_undir_br02, ne_inf_undir_br02, 'o', label='BR02', markersize=14)
plt.plot(ne_sal_undir_br26, ne_inf_undir_br26, 'o', label='BR26', markersize=14)
plt.plot(ne_sal_undir_or13, ne_inf_undir_or13, 'o', label='OR13', markersize=14)
plt.plot(ne_sal_undir_or02, ne_inf_undir_or02, 'o', label='OR02', markersize=14)
plt.plot(ne_sal_undir_wh27, ne_inf_undir_wh27, 'o', label='WH27', markersize=14, color='m')
plt.plot(ne_sal_undir_wh57, ne_inf_undir_wh57, 'o', label='WH57', markersize=14, color='c')
plt.plot(ne_sal_undir_y437, ne_inf_undir_y437, 'o', label='Y437', markersize=14)
plt.plot(ne_sal_undir_br0, ne_inf_undir_br0, 'o', label='BR0', markersize=14)
plt.plot([0.4, 0.9], [0.4, 0.9], 'k--')
plt.xlim([0.5, 0.9])
plt.ylim([0.5, 0.9])
plt.legend(fontsize=14)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.title('Spectral Variability (A.U.) - Undirected Song', fontsize=20)
plt.ylabel('Norepinephrine',fontsize=20)
plt.xlabel('Saline',fontsize=20)
plt.savefig('sv_undirected_ne.png', dpi=300)

# Directed NE Infusion vs. Saline
plt.figure(figsize=(12,8))
plt.plot(ne_sal_dir_br02, ne_inf_dir_br02, 'o', label='BR02', markersize=14)
plt.plot(ne_sal_dir_br26, ne_inf_dir_br26, 'o', label='BR26', markersize=14)
plt.plot(ne_sal_dir_or13, ne_inf_dir_or13, 'o', label='OR13', markersize=14)
plt.plot(ne_sal_dir_or02, ne_inf_dir_or02, 'o', label='OR02', markersize=14)
plt.plot(ne_sal_dir_wh27, ne_inf_dir_wh27, 'o', label='WH27', markersize=14, color = 'm')
plt.plot(ne_sal_dir_wh57, ne_inf_dir_wh57, 'o', label='WH57', markersize=14, color = 'c')
plt.plot([0.4, 0.9], [0.4, 0.9], 'k--')
plt.xlim([0.5, 0.9])
plt.ylim([0.5, 0.9])
plt.legend(fontsize=14)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.title('Spectral Variability (A.U.) - Directed Song', fontsize=20)
plt.ylabel('Norepinephrine', fontsize=20)
plt.xlabel('Saline', fontsize=20)
plt.savefig('sv_directed_ne.png', dpi=300)

# PHE song frequency

# Directed PHE Infusion vs. Saline
plt.figure(figsize=(12,8))
plt.plot(phe_sal_dir_or46, phe_inf_dir_or46, 'o', label='OR46', markersize=14)
plt.plot(phe_sal_dir_br26, phe_inf_dir_br26, 'o', label='BR26', markersize=14)
plt.plot(phe_sal_dir_wh27, phe_inf_dir_wh27, 'o', label='WH27', markersize=14, color='m')
plt.plot(phe_sal_dir_wh57, phe_inf_dir_wh57, 'o', label='WH57', markersize=14, color='c')
plt.plot([0.4, 0.9], [0.4, 0.9], 'k--')
plt.xlim([0.5, 0.9])
plt.ylim([0.5, 0.9])
plt.legend(fontsize=14)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.title('Spectral Variability (A.U.) - Directed Song', fontsize=20)
plt.ylabel('Phentolamine', fontsize=20)
plt.xlabel('Saline ', fontsize=20)
plt.savefig('sv_directed_phe.png', dpi=300)

# Box and whisker plot - LC Stimulation
lc_rd08_overall = [lc_stim_dir_rd08, lc_nostim_dir_rd08, lc_stim_undir_rd08, lc_nostim_undir_rd08]
lc_wh09_overall = [lc_stim_dir_wh09, lc_nostim_dir_wh09, lc_stim_undir_wh09, lc_nostim_undir_wh09]

# LC boxplot
lc_overall = [lc_stim_undir_overall, lc_nostim_undir_overall]
plt.figure(figsize=(12,8))
ax1 = sns.boxplot(data=lc_overall, showfliers=False, color='c')
ax2 = sns.swarmplot(data=lc_overall, color=".25", size=14)
ax1.set_xticklabels(['LC Stim', 'No Stim'], fontsize=20)
plt.ylabel('Spectral Variability (A.U.)', fontsize=20)
plt.title('Undirected Song', fontsize=20)
plt.ylim([0.4, 0.9])
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.tight_layout()

lc_undir_x1, lc_undir_x2 = 0, 1
lc_undir_y, h, col = 0.865, 0.01, 'k'
plt.plot([lc_undir_x1, lc_undir_x1, lc_undir_x2, lc_undir_x2], [lc_undir_y, lc_undir_y+h, lc_undir_y+h, lc_undir_y], lw=1.5, c=col)
plt.text((lc_undir_x1+lc_undir_x2)*.5, lc_undir_y+h, r'* ($p = 3.5 \cdot 10^{-6}$)', ha='center', va='bottom', color=col, fontsize=18);
plt.savefig('boxwhisker_undirected_lc.png', dpi=300)

# box and whisker plot - NE Infusion
ne_overall = [ne_inf_dir_overall, ne_sal_dir_overall, ne_inf_undir_overall, ne_sal_undir_overall]
plt.figure(figsize=(12,8))
ax1 = sns.boxplot(data=ne_overall[2:], showfliers=False, color='c')
ax2 = sns.swarmplot(data=ne_overall[2:], color=".25", size=14)
plt.title('Undirected Song', fontsize=20)
plt.ylabel('Spectral Variability (A.U.)', fontsize=20)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.ylim([0.4, 1.0])
ax1.set_xticklabels(['Norepinephrine', 'Saline'], fontsize=20)
plt.tight_layout()
undir_x1, undir_x2 = 0, 1
undir_y, h, col = 0.8, 0.025, 'k'
plt.plot([undir_x1, undir_x1, undir_x2, undir_x2], [undir_y, undir_y+h, undir_y+h, undir_y], lw=1.5, c=col)
plt.text((undir_x1+undir_x2)*.5, undir_y+h, r'* ($p = 8.0 \cdot 10^{-7}$)', ha='center', va='bottom', color=col, fontsize=18);
plt.savefig('boxwhisker_ne.png', dpi=300)

# phe boxplot
phe_overall = [phe_inf_overall, phe_sal_overall]
plt.figure(figsize=(12,8))
ax1 = sns.boxplot(data=phe_overall, showfliers=False, color='c')
ax2 = sns.swarmplot(data=phe_overall, color=".25", size=14)
plt.title('Directed Song', fontsize=20)
plt.ylabel('Spectral Variability (A.U.)', fontsize=20)
plt.ylim([0.4, 1.0])
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
ax1.set_xticklabels(['Phentolamine', 'Saline'], fontsize=20)
plt.tight_layout()

phe_dir_x1, phe_dir_x2 = 0, 1
phe_dir_y, h, col = 0.9, 0.025, 'k'
plt.plot([phe_dir_x1, phe_dir_x1, phe_dir_x2, phe_dir_x2], [phe_dir_y, phe_dir_y+h, phe_dir_y+h, phe_dir_y], lw=1.5, c=col)
plt.text((phe_dir_x1+phe_dir_x2)*.5, phe_dir_y+h, r'* ($p = 2.6 \cdot 10^{-4}$)', ha='center', va='bottom', color=col, fontsize=18);

plt.savefig('phe_directed_boxplot.png', dpi=300)

# song frequency plots

plt.figure(figsize=(12,8))
plt.plot(rd02_songfreq_undir_nostim, rd02_songfreq_undir_lcstim, 'o', markersize=14, label='RD02')
plt.plot(pu31_songfreq_undir_nostim, pu31_songfreq_undir_lcstim, 'o', markersize=14, label='PU31')
plt.plot(bl16_songfreq_undir_nostim, bl16_songfreq_undir_lcstim, 'o', markersize=14, label='BL16')
plt.plot(si026_songfreq_undir_nostim, si026_songfreq_undir_lcstim, 'o', markersize=14, label='SI026')
plt.plot(rd08_songfreq_undir_nostim, rd08_songfreq_undir_lcstim, 'o', markersize=14, label='RD08')
plt.plot(wh09_songfreq_undir_nostim, wh09_songfreq_undir_lcstim, 'ko', fillstyle='none', markersize=10, label='WH09')
plt.plot([0, 3], [0, 3], 'k--')
plt.title('Singing Frequency (songs/min) - Undirected Song', fontsize=20)
plt.xlabel('No Stim', fontsize=20)
plt.xlim([0, 2.5])
plt.ylabel('LC Stim', fontsize=20)
plt.ylim([0, 2.5])
plt.tight_layout()
plt.legend(fontsize=14)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.savefig('songfreq_lc_zach.png', dpi=300)

plt.figure(figsize=(12,8))
plt.plot(br02_songfreq_undir_sal, br02_songfreq_undir_ne, 'o', markersize=14, label='BR02')
plt.plot(br26_songfreq_undir_sal, br26_songfreq_undir_ne, 'o', markersize=14, label='BR26')
plt.plot(or13_songfreq_undir_sal, or13_songfreq_undir_ne, 'o', markersize=14, label='OR13')
plt.plot(or02_songfreq_undir_sal, or02_songfreq_undir_ne, 'o', markersize=14, label='OR02')
plt.plot(wh27_songfreq_undir_sal, wh27_songfreq_undir_ne, 'o', markersize=14, label='WH27')
plt.plot(wh57_songfreq_undir_sal, wh57_songfreq_undir_ne, 'o', markersize=14, label='WH57')
plt.plot(y437_songfreq_undir_sal, y437_songfreq_undir_ne, 'o', markersize=14, label='Y437')
plt.plot(br0_songfreq_undir_sal, br0_songfreq_undir_ne, 'o', markersize=14, label='BR0')
plt.plot([0, 4], [0, 4], 'k--')
plt.title('Singing Frequency (songs/min) - Undirected Song', fontsize=20)
plt.xlabel('Saline', fontsize=20)
plt.xlim([0, 4])
plt.ylabel('Norepinephrine', fontsize=20)
plt.ylim([0, 4])
plt.tight_layout()
plt.legend(fontsize=14)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.savefig('songfreq_ne_zach.png', dpi=300)

plt.figure(figsize=(12,8))
plt.plot(or46_songfreq_dir_sal, or46_songfreq_dir_phe, 'o', markersize=14, label='OR46')
plt.plot(br26_songfreq_dir_sal, br26_songfreq_dir_phe, 'o', markersize=14, label='BR26')
plt.plot(wh27_songfreq_dir_sal, wh27_songfreq_dir_phe, 'o', markersize=14, label='WH27')
plt.plot(wh57_songfreq_dir_sal, wh57_songfreq_dir_phe, 'o', markersize=14, label='WH57')
plt.plot([0, 7], [0, 7], 'k--')
plt.title('Singing Frequency (songs/min) - Directed Song', fontsize=20)
plt.xlabel('Saline', fontsize=20)
plt.xlim([0, 7])
plt.ylabel('Phentolamine', fontsize=20)
plt.ylim([0, 7])
plt.tight_layout()
plt.legend(fontsize=14)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.savefig('songfreq_phe_zach.png', dpi=300)

ne_averages = []
for i in range(0, len(ne_overall)):
    ne_averages.append(np.sum(ne_overall[i])/len(ne_overall[i]))

ne_dir_t_val, ne_dir_p_val = ttest_rel(ne_overall[0], ne_overall[1])
ne_dir_p_val = ne_dir_p_val / 2.0
ne_undir_t_val, ne_undir_p_val = ttest_rel(ne_overall[2], ne_overall[3])
ne_undir_p_val = ne_undir_p_val / 2.0
ne_inf_t_val, ne_inf_p_val = ttest_rel(ne_overall[0], ne_overall[2][:len(ne_overall[0])])
ne_inf_p_val = ne_inf_p_val / 2.0
ne_sal_t_val, ne_sal_p_val = ttest_rel(ne_overall[1], ne_overall[3][:len(ne_overall[1])])
ne_sal_p_val = ne_sal_p_val / 2.0

rd08_dir_t_val, rd08_dir_p_val_twoTailed = ttest_rel(lc_rd08_overall[0], lc_rd08_overall[1])
rd08_dir_p_val_oneTailed = rd08_dir_p_val_twoTailed / 2.0
rd08_undir_t_val, rd08_undir_p_val = ttest_rel(lc_rd08_overall[2], lc_rd08_overall[3])
rd08_undir_p_val = rd08_undir_p_val / 2.0

lc_undir_t_val, lc_undir_p_val = ttest_rel(lc_overall[0], lc_overall[1])
lc_undir_p_val = lc_undir_p_val / 2.0

phe_dir_t_val, phe_dir_p_val = ttest_rel(phe_inf_overall, phe_sal_overall)
phe_dir_p_val= phe_dir_p_val / 2.0

wh09_dir_t_val, wh09_dir_p_val = ttest_rel(lc_wh09_overall[0], lc_wh09_overall[1])
wh09_dir_p_val = wh09_dir_p_val / 2.0
wh09_undir_t_val, wh09_undir_p_val = ttest_rel(lc_stim_undir_wh09, lc_nostim_undir_wh09)
wh09_undir_p_val = wh09_undir_p_val / 2.0

# some old(?) data I found in Chris' excel files for PHE directed
dir_phe_updated=[0.671,0.668,0.627,0.745,0.754,0.735,0.645,0.778,0.664,0.881,0.862,0.885,0.890,0.916,0.840,0.921,0.760,0.729,0.825,0.867,0.903,0.915]
dir_sal_updated = [0.647,0.658,0.594,0.716,0.730,0.728,0.658,0.764,0.630,0.866,0.838,0.863,0.852,0.900,0.841,0.916,0.777,0.782,0.881,0.857,0.965,0.916]

updated_phe_dir_t_val, updated_phe_p_val = ttest_rel(dir_phe_updated, dir_sal_updated)
updated_phe_p_val = updated_phe_p_val/ 2.0

######## boxplots for individual birds ###########
import matplotlib
matplotlib.rcParams.update({'font.size': 18})

# plot individual syllables for each bird and make separate plots for each of them
def plot_individual_sylls(control,exp,birdnames,directed,conds_str,save_filename):
    if directed:
        title = 'Directed Song'
    else:
        title = 'Undirected Song'
    for idx,(c,e) in enumerate(zip(control, exp)):
        fig,ax = plt.subplots(1,1,figsize=(10,10))
        # plot individual syllable lines between no stim and stim
        for y_c,y_e in zip(c,e):
            plt.plot([0.02,0.98], [y_c,y_e],color='k',linestyle='--')
        # plot data points
        plt.scatter(np.zeros(len(c)),c,color='blue',s=[100 for pt in c])
        plt.scatter(np.ones(len(e)),e,color='red',s=[100 for pt in e])
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.xaxis.set_ticks_position('none')
        ax.yaxis.set_ticks_position('none')
        plt.xlim([-0.25, 1.25])
        plt.xticks([0,1],tuple(conds_str))
        plt.ylabel('Spectral Variability (A. U.)')
        plt.ylim([0.4,0.9])
        plt.title('{} - {}'.format(title,birdnames[idx]))
        plt.savefig(save_filename+'_{}.png'.format(birdnames[idx]),dpi=300)

# plot individual syllables for each bird and amalgamate into one big figure
def plot_overall_sylls(control,exp,birdnames,directed,conds_str,save_filename):
    if directed:
        title = 'Directed Song'
    else:
        title = 'Undirected Song'

    x_indices_left = np.arange(0,len(birdnames)*2,2)

    fig,ax = plt.subplots(1,1,figsize=(18,10))
    for (c,e,idx_left) in zip(control, exp,x_indices_left):
        # only include label for one instance so legend isn't redundant
        if idx_left == x_indices_left[0]:
            plt.scatter(np.random.normal(idx_left,0,size=len(c)),c,color='blue',label=conds_str[0],s=[120 for pt in c])
            plt.scatter(np.random.normal(idx_left+1,0,size=len(e)),e,color='red',label=conds_str[1],s=[120 for pt in e])
        else:
            plt.scatter(np.random.normal(idx_left,0,size=len(c)),c,color='blue',s=[120 for pt in c])
            plt.scatter(np.random.normal(idx_left+1,0,size=len(e)),e,color='red',s=[120 for pt in e])
        # plot individual syllable lines between no stim and stim
        for y_c,y_e in zip(c,e):
            plt.plot([idx_left+0.04,idx_left+0.96], [y_c,y_e],color='k',linestyle='--')
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.xaxis.set_ticks_position('none')
    ax.yaxis.set_ticks_position('none')
    plt.xticks([x+0.5 for x in x_indices_left],tuple(birdnames))
    plt.ylabel('Spectral Variability (A. U.)')
    plt.ylim([0.4,0.9])
    plt.title(title)
    plt.legend()
    plt.savefig(save_filename+'.png',dpi=300)

# Undirected - LC Stim vs. No Stim
lc_nostim_undir = [lc_nostim_undir_rd02, lc_nostim_undir_pu31, lc_nostim_undir_bl16, lc_nostim_undir_si26,lc_nostim_undir_wh09]
lc_stim_undir = [lc_stim_undir_rd02, lc_stim_undir_pu31, lc_stim_undir_bl16, lc_stim_undir_si26,lc_stim_undir_wh09]
birdnames = ['RD02', 'PU31', 'BL16', 'SI26', 'WH09 (Control)']

#plot_individual_sylls(lc_nostim_undir,lc_stim_undir,birdnames,False,['No Stim', 'Stim'],'lc_undir_syllables')
plot_overall_sylls(lc_nostim_undir,lc_stim_undir,birdnames,False,['No Stim', 'Stim'],'lc_undir_syllables')


# Undirected - NE vs. Saline
ne_inf_undir = [ne_inf_undir_br02,ne_inf_undir_br0,ne_inf_undir_or13,ne_inf_undir_br26,ne_inf_undir_or02,ne_inf_undir_wh27,ne_inf_undir_wh57,ne_inf_undir_y437]
ne_sal_undir = [ne_sal_undir_br02,ne_sal_undir_br0,ne_sal_undir_or13,ne_sal_undir_br26,ne_sal_undir_or02,ne_sal_undir_wh27,ne_sal_undir_wh57,ne_sal_undir_y437]
birdnames_ne_undir = ['BR02', 'BR00', 'OR13', 'BR26', 'OR02', 'WH27', 'WH57', 'YW437']

#plot_individual_sylls(ne_sal_undir,ne_inf_undir,birdnames_ne_undir,False,['Saline', 'NE'],'ne_undir_syllables')
plot_overall_sylls(ne_sal_undir,ne_inf_undir,birdnames_ne_undir,False,['Saline', 'NE'],'ne_undir_syllables')

# Directed - NE vs. Saline
ne_inf_dir = [ne_inf_dir_br02,ne_inf_dir_or13,ne_inf_dir_br26,ne_inf_dir_or02,ne_inf_dir_wh27,ne_inf_dir_wh57]
ne_sal_dir = [ne_sal_dir_br02,ne_sal_dir_or13,ne_sal_dir_br26,ne_sal_dir_or02,ne_sal_dir_wh27,ne_sal_dir_wh57]
birdnames_ne_dir = ['BR02', 'OR13', 'BR26', 'OR02', 'WH27', 'WH57']

#plot_individual_sylls(ne_sal_dir,ne_inf_dir,birdnames_ne_dir,True,['Saline', 'NE'],'ne_dir_syllables')
plot_overall_sylls(ne_sal_dir,ne_inf_dir,birdnames_ne_dir,True,['Saline', 'NE'],'ne_dir_syllables')

# Directed - PHE vs. Saline
phe_inf_dir = [phe_inf_dir_or46,phe_inf_dir_br26,phe_inf_dir_wh27,phe_inf_dir_wh57]
phe_sal_dir = [phe_sal_dir_or46,phe_sal_dir_br26,phe_sal_dir_wh27,phe_sal_dir_wh57]
birdnames_phe_dir = ['OR46', 'BR26', 'WH27', 'WH57']

#plot_individual_sylls(phe_sal_dir,phe_inf_dir,birdnames_phe_dir,True,['Saline', 'PHE'],'phe_dir_syllables')
plot_overall_sylls(phe_sal_dir,phe_inf_dir,birdnames_phe_dir,True,['Saline', 'PHE'],'phe_dir_syllables')
