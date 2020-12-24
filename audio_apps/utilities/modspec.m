function z = modspec(x,y)

z = mod(x,y);
z(z==0) = y; 