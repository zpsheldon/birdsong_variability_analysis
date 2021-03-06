function h = plotDTW2(map,w,x,y,dt)

h = zeros(1,3);

tx = 1:size(x,2);
ty = 1:size(y,2);
f = 1:size(x,1);

h(1) = subplot(4,4,[2:4,6:8,10:12]);

clims = [0 prctile(map(:),95)];

imagesc(tx,ty,map',clims);
set(gca,'ydir','normal');
hold on
plot(w(:,1),w(:,2),'w','linewidth',1.5);
set(gca,'ylim',[min(ty) max(ty)],'xlim',[min(tx) max(tx)])



set(gca,'xtick',[],'ytick',[]);

if size(x,1)>1
    h(2) = subplot(4,4,[1 5 9]);
    imagesc(f,ty,y');
    set(gca,'ylim',[min(ty) max(ty)],'ydir','normal')
    
    
    h(3) = subplot(4,4,14:16);
    
    imagesc(tx,f,x);
    set(gca,'ydir','normal');
    set(gca,'xlim',[min(tx) max(tx)])
    
else
    
    h(2) = subplot(4,4,[1 5 9]);
    plot(y,ty,'k-','linewidth',2);
    set(gca,'ylim',[min(ty) max(ty)],'ydir','normal')
    
    
    h(3) = subplot(4,4,14:16);
    
    plot(tx,x','k-','linewidth',2);
    set(gca,'ydir','normal');
    set(gca,'xlim',[min(tx) max(tx)])
    
    
end

