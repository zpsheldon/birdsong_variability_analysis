function Y = chunk(X,v)


chunknm = length(v)-1;
Y = zeros(size(X,1),chunknm);

for i = 1:chunknm
    
    if i < chunknm
        indstmp = v(i):v(i+1)-1;
    else
        indstmp = v(i):v(i+1);
    end

    if size(indstmp,2) > 1
        Y(:,i) = sum(X(:,indstmp)')';
    else
        Y(:,i) = X(:,indstmp);
    end
end