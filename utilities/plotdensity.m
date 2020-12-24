function h = plotdensity(X,binm)

if nargin == 1
    binm = 50;
end

[n,d] = size(X);

h = figure;

subind = 1;

for d1 = 1:d
    for d2 = d1+1:d
        
        subplot(floor(d/2),ceil((d*(d-1)/2)/floor(d/2)),subind)
        
        lim11 = prctile(X(:,d1),0);
        lim12 = prctile(X(:,d1),100);
        
        dlim = (prctile(X(:,d1),75) - prctile(X(:,d1),25))/binm;

        vc1 = lim11:dlim:lim12;
        
        lim21 = prctile(X(:,d2),0);
        lim22 = prctile(X(:,d2),100);
        
        
        dlim = (prctile(X(:,d2),75) - prctile(X(:,d2),25))/binm;
        
        vc2 = lim21:dlim:lim22;
        
        
        [N,C] = hist3(X(:,[d1 d2]),{vc1,vc2});
        imagesc(C{1},C{2},N');
        
%         xlim([prctile(X(:,d1),1) prctile(X(:,d1),99)])
%         ylim([prctile(X(:,d2),1) prctile(X(:,d2),99)])
%          
        hold on
        
        set(gca,'ydir','normal')
        
        subind = subind + 1;
        
    end
end
        