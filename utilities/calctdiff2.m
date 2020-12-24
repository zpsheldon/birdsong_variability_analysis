function overlap = calctdiff2(tmvc,thresh)

overlap = 0;

for i = 1:length(tmvc)
    for j = i+1:length(tmvc)
        
        if abs(tmvc(i)-tmvc(j))<=thresh
           overlap = overlap + 1; 
        end
        
    end
end