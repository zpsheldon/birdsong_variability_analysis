function [l,u] = propCI(p,n,alpha)

q = 1 - p;
z = norminv(1-alpha/2,0,1);

l = (p + z^2 ./ (2*n) - z * sqrt(p.*q./n + (z^2) ./(4*n.^2))) ./ (1 + z^2 ./ n);
u = (p + z^2 ./ (2*n) + z * sqrt(p.*q./n + (z^2) ./(4*n.^2))) ./ (1 + z^2 ./ n);