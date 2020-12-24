function plotgauss(h,mu,sigma)

x1 = -1:.01:1;
y1 = sqrt(1-x1.^2);
x1 = [x1 fliplr(x1)];
y1 = [y1 fliplr(-y1)];

U = [x1' y1'];

k = size(mu,1);

for classind = 1:k
    
    sigmatmp = squeeze(sigma(:,:,classind));
    mutmp = mu(classind,:);
    
    Schol = chol(sigmatmp);
    W = 1.96*U*Schol;
    
    plot(h,W(:,1)+mutmp(1),W(:,2)+mutmp(2),'r-','linewidth',2);
    
end