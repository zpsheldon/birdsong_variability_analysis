function overlap = calctdiff(tms1,tms2,thresh)

len1 = length(tms1);
len2 = length(tms2);

overlap = 0;

for i = 1:len1
    for j = 1:len2
        
        if abs(tms1(i)-tms2(j))<=thresh
           overlap = overlap + 1; 
        end
        
    end
end