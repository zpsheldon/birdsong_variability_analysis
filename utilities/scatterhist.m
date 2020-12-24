function [h,tnorm,lab] = scatterhist(x,y,horg,bin)

if nargin < 4
    bin = 10;
end

dx = range(x)/bin;
dy = range(y)/bin;

x = round(x/dx)*dx;
y = round(y/dy)*dy;

[t,p1,p2,lab] = crosstab(x,y);
xvals = sort(unique(x));
yvals = sort(unique(y));

xmat = repmat(xvals',[1 length(yvals)]);
ymat = repmat(yvals,[length(xvals) 1]);

t = max(t,.0001); %*range(xvals)*range(yvals);

if nargin == 3
    h = scatter(horg,xmat(:),ymat(:),sqrt(t(:)));
else
    h = scatter(xmat(:),ymat(:),t(:));
end
