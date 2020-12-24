function d = Afit(A,dt,coeffs,const)

m = size(A,1);
d = zeros(m);
n = size(coeffs,1);

for i = 1:n
   d = d + (A^dt(i))*squeeze(coeffs(i,:,:));
end

d = sum(sum((d-const).^2));