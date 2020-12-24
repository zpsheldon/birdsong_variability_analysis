function dhatdist = bootstrp_grouped_dist(d,grp,bootN,func_handle)

grpU = unique(grp);
grpN = numel(grpU);

T = tabulate(grp);
grpU = T(:,1);
grpM = cell2mat(T(:,2));

dhatdist = zeros(bootN,1);

for bootind = 1:bootN
    grpinds = ceil(rand(grpN,1)*grpN);
    
    Ntmp = sum(grpM(grpinds));
    
    k = 1;
    dtmp = zeros(Ntmp,1);
    
    for grpind = 1:numel(grpinds)   
        indstmp = ismember(grp,grpU(grpinds(grpind)));
        indstmp2 = k:k+sum(indstmp)-1;
        
        bootK = sum(indstmp);
        bootindstmp = ceil(bootK*rand(bootK,1));
        bootindstmp2 = find(indstmp);
        bootindstmp2 = bootindstmp2(bootindstmp);
        
        dtmp(indstmp2) = d(bootindstmp2);
        
        k = k + sum(indstmp);
    end
    
    dhatdist(bootind) = func_handle(dtmp);

end

