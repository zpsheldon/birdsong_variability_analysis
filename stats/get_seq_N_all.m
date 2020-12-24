clear
clc
close all

birdarrDIR = {'BR_2','OR_13','BR26','OR2','WH27','OR_46'};
birdarrLC = {'RD_2','SI_026','PU_31','BL_16'};
birdarrPHE = {'OR_46','BR26','WH27','WH57'};
birdarrNE =  {'BR_2','OR_13','BR26','OR2','WH57','Y437','WH27','BRFINAL'};

paramatNE = get_seq_N(birdarrNE,1);
paramatDIR = get_seq_N(birdarrDIR,2);
paramatLC = get_seq_N(birdarrLC,0);
paramatPHE = get_seq_N(birdarrPHE,3);
    
    

paramat_all = [paramatNE;paramatDIR;paramatLC;paramatPHE]