function h = plotDTW(h_parent,map,w,x,y,t,f)

h = zeros(1,3);

tx = 1:size(x,2);
ty = 1:size(y,2);
f = 1:size(x,1);

h(1) = subplot(4,4,[2:4,6:8,10:12],'Parent',h_parent);

imagesc(tx,ty,map');
set(gca,'ydir','normal');
hold on
plot(w(:,1),w(:,2),'k','linewidth',1.5);
set(gca,'ylim',[min(ty) max(ty)],'xlim',[min(tx) max(tx)])


set(gca,'xtick',[],'ytick',[]);


h(2) = subplot(4,4,[1 5 9],'Parent',h_parent);
imagesc(f,ty,y');
set(gca,'ylim',[min(ty) max(ty)],'ydir','normal')


h(3) = subplot(4,4,14:16,'Parent',h_parent);

imagesc(tx,f,x);
set(gca,'ydir','normal');
set(gca,'xlim',[min(tx) max(tx)])

