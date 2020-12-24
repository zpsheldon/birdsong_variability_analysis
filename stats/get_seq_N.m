function paramat = get_seq_N(birdarr,cond)

Nsal = zeros(numel(birdarr),1);
Nexp = Nsal;
Nsyl = Nsal;

for birdind = 1:numel(birdarr)
    
    birdind
 
%     load(['./' birdarr{birdind} '/timedata.mat'])
%     
%     switch cond
%         case 1
%             lenseqexp = timestrct_U_NE.lenseq;
%             lenseqsal = timestrct_U_sal.lenseq;
%             seq = timestrct_U_sal.seqarr;
%             
%         case 0
%             lenseqexp = timestrct_U_stim.lenseq;
%             lenseqsal = timestrct_U_nostim.lenseq;
%             seq = timestrct_U_nostim.seqarr;
%             
%         case 2
%             lenseqexp = timestrct_F_sal.lenseq;
%             lenseqsal = timestrct_U_sal.lenseq;
%             seq = timestrct_U_sal.seqarr;
%             
%         case 3
%             lenseqsal = timestrct_F_sal.lenseq;
%             lenseqexp = timestrct_F_PHE.lenseq;
%             seq = timestrct_F_sal.seqarr;
%             
%     end
    

    load(['./' birdarr{birdind} '/specdata.mat'])
    
    switch cond
        case 1
            lenseqexp = specstrct_U_NE.distmat;
            lenseqsal = specstrct_U_sal.distmat;
            
        case 0
            lenseqexp = specstrct_U_stim.distmat;
            lenseqsal = specstrct_U_nostim.distmat;
            
        case 2
            lenseqexp = specstrct_F_sal.distmat;
            lenseqsal = specstrct_U_sal.distmat;
            
        case 3
            lenseqsal = specstrct_F_sal.distmat;
            lenseqexp = specstrct_F_PHE.distmat;

            
    end


    Nsal(birdind) = size(lenseqsal,1);
    Nexp(birdind) = size(lenseqexp,1);

    Nsyl(birdind) = size(lenseqsal,2);

    
end


paramat = [Nsal Nexp Nsyl];