function Y = chunk2(X,sg)


chunknm = sum(diff(sg)~=0)+1;
Y = zeros(size(X,1),chunknm);

k = 1;

for i = 1:chunknm    
    if k<=length(sg) && sg(k)==1
        while k<=length(sg) && sg(k)==1
            Y(:,i) = Y(:,i)+X(:,k);
            k = k + 1;
        end
    else
        Y(:,i) = X(:,k);
        k = k + 1;
    end
end