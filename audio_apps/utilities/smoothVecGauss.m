function smoothVec = smoothVecGauss(vec,n)

b = gausswin(n)';
b = b / sum(b);

[rowNum,colNum] = size(vec);

smoothVec = zeros(rowNum,colNum-n+1);

for row = 1:rowNum
    smoothVec(row,:) = xcorr212(b,vec(row,:),0);
%     for col = 1:colNum - n + 1
%         smoothVec(row,col) = sum(vec(row,col:col + n - 1).*b);
%     end
end