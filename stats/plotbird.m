function h = plotbird(x,y,plotbirdarr,ax)

load birdplotparams.mat birdarr markerarr colorarr sizearr


hold(ax,'on')
xN = numel(x);

for xind = 1:xN
    
    birdind = find(ismember(birdarr,plotbirdarr{xind}));
    
    marker = markerarr{birdind};
    clr = colorarr{birdind};
    sz = sizearr(birdind)*.75;
    
    plot(x(xind),y(xind),'marker',marker,'color',clr,'markersize',sz,'linewidth',2)
    
end

set(gca,'box','off','ticklength',[0,0],'fontsize',14)

h = gca