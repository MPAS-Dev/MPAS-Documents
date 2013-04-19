function [ area ] = areaBasedonVertexCoords( tri )
%areaBasedonVertexCoords : given vertex coord, compute area
%coords need to be listed in CCW direction to get positive area

x1=tri(1,1);
x2=tri(2,1);
x3=tri(3,1);

y1=tri(1,2);
y2=tri(2,2);
y3=tri(3,2);

area = 0.5*((x1-x2)*(y1-y3)-(y1-y2)*(x1-x3));


end

