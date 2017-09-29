function [ u, v ] = testVector( x,y )
%testVector: given x and y, produce a test vector

u0 = 1.0;
ux = 0.0;
uy = -1.0;

v0 = 1.0;
vx = 1.0;
vy = 0.0;

u = u0 + ux*x + uy*y;
v = v0 + vx*x + vy*y;

end

