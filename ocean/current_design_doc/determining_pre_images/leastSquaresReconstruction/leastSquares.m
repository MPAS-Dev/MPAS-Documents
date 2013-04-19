clear all

%testing the least squares approach to vector reconstruction

nVertices = 10;
xVertex = zeros(nVertices,1);
yVertex = zeros(nVertices,1);

nEdges = 9;
xEdge = zeros(nEdges,1);
yEdge = zeros(nEdges,1);
uVector = zeros(nEdges,2);
normal = uVector;

% everything between the ---- lines just builds the vertices and edges
%{x,y}Vertex: position of vertices
%{x,y}Edge: position of edges
%normal: unit vector normal to edges
%uVector: test vector evaluated at (xEdge,yEdge)
%------------------------------------------------------------------------
dx = 1.0;
dtr = pi/180.0;
delta = 120.0*dtr;

xVertex(1)=0.0;
yVertex(1)=0.0;

angle = 0.0*dtr;
iVertexBase=1;
iCount = 1;
iEdge = 0;
for i=1:3
    iCount=iCount+1;
    xVertex(iCount) = xVertex(iVertexBase)+dx*cos(angle);
    yVertex(iCount) = yVertex(iVertexBase)+dx*sin(angle);
    iEdge=iEdge+1;
    xEdge(iEdge) = 0.5*(xVertex(iVertexBase)+xVertex(iCount));
    yEdge(iEdge) = 0.5*(yVertex(iVertexBase)+yVertex(iCount));
    [uVector(iEdge,1), uVector(iEdge,2)] = ...
                          testVector(xEdge(iEdge),yEdge(iEdge));
    normal(iEdge,2) = +(xVertex(iCount) - xVertex(iVertexBase));
    normal(iEdge,1) = -(yVertex(iCount) - yVertex(iVertexBase));
    angle = angle + delta;
end

angle = 300.0*dtr;
iVertexBase=2;
for i=1:2
    iCount=iCount+1;
    xVertex(iCount) = xVertex(iVertexBase)+dx*cos(angle);
    yVertex(iCount) = yVertex(iVertexBase)+dx*sin(angle);
    iEdge=iEdge+1;
    xEdge(iEdge) = 0.5*(xVertex(iVertexBase)+xVertex(iCount));
    yEdge(iEdge) = 0.5*(yVertex(iVertexBase)+yVertex(iCount));
    [uVector(iEdge,1), uVector(iEdge,2)] = ...
                          testVector(xEdge(iEdge),yEdge(iEdge));
    normal(iEdge,2) = +(xVertex(iCount) - xVertex(iVertexBase));
    normal(iEdge,1) = -(yVertex(iCount) - yVertex(iVertexBase));
    angle = angle + delta;
end

angle = 60.0*dtr;
iVertexBase=3;
for i=1:2
    iCount=iCount+1;
    xVertex(iCount) = xVertex(iVertexBase)+dx*cos(angle);
    yVertex(iCount) = yVertex(iVertexBase)+dx*sin(angle);
    iEdge=iEdge+1;
    xEdge(iEdge) = 0.5*(xVertex(iVertexBase)+xVertex(iCount));
    yEdge(iEdge) = 0.5*(yVertex(iVertexBase)+yVertex(iCount));
    [uVector(iEdge,1), uVector(iEdge,2)] = ...
                          testVector(xEdge(iEdge),yEdge(iEdge));
    normal(iEdge,2) = +(xVertex(iCount) - xVertex(iVertexBase));
    normal(iEdge,1) = -(yVertex(iCount) - yVertex(iVertexBase));
    angle = angle + delta;
end

angle = 180.0*dtr;
iVertexBase=4;
for i=1:2
    iCount=iCount+1;
    xVertex(iCount) = xVertex(iVertexBase)+dx*cos(angle);
    yVertex(iCount) = yVertex(iVertexBase)+dx*sin(angle);
    iEdge=iEdge+1;
    xEdge(iEdge) = 0.5*(xVertex(iVertexBase)+xVertex(iCount));
    yEdge(iEdge) = 0.5*(yVertex(iVertexBase)+yVertex(iCount));
    [uVector(iEdge,1), uVector(iEdge,2)] = ...
                          testVector(xEdge(iEdge),yEdge(iEdge));
    normal(iEdge,2) = +(xVertex(iCount) - xVertex(iVertexBase));
    normal(iEdge,1) = -(yVertex(iCount) - yVertex(iVertexBase));
    angle = angle + delta;
end

for iEdge=1:nEdges
    mag = sqrt( normal(iEdge,1)^2 + normal(iEdge,2)^2);
    normal(iEdge,:) = normal(iEdge,:) / mag;
end
%------------------------------------------------------------------------

%set aside space for matrix and rhs
M = zeros(nEdges,6);
rhs = zeros(nEdges,1);

%fill matrix and rhs
%rhs is just the test vector dotted with the normal
%Matrix(:,1) is n_x
%Matrix(:,2) is n_x*x
%Matrix(:,3) is n_x*y
%Matrix(:,4) is n_x
%Matrix(:,5) is n_y*x
%Matrix(:,6) is n_y*y

for iEdge=1:nEdges
    rhs(iEdge) = normal(iEdge,1)*uVector(iEdge,1) + ...
                 normal(iEdge,2)*uVector(iEdge,2);

    M(iEdge,1) = normal(iEdge,1);
    M(iEdge,2) = normal(iEdge,1)*xEdge(iEdge);
    M(iEdge,3) = normal(iEdge,1)*yEdge(iEdge);
    
    M(iEdge,4) = normal(iEdge,2);
    M(iEdge,5) = normal(iEdge,2)*xEdge(iEdge);
    M(iEdge,6) = normal(iEdge,2)*yEdge(iEdge);
    
end

%find psuedo inversion
M = pinv(M);

%solve!
solution = M*rhs

%put solution and uVector into single array so they are plotted to scale
nTmp = nEdges+1;
xTmp = zeros(nTmp,1);
yTmp = zeros(nTmp,1);
uTmp = zeros(nTmp,2);
uTmp(1:nEdges,:) = uVector(:,:);
xTmp(1:nEdges) = xEdge(:);
yTmp(1:nEdges) = yEdge(:);
uTmp(nTmp,1) = solution(1);
uTmp(nTmp,2) = solution(4);

scatter(xVertex,yVertex)
axis([-3 3 -3 3])
axis equal
hold on
scatter(xEdge,yEdge,'filled')
quiver(xTmp,yTmp, uTmp(:,1),uTmp(:,2))
%quiver(xEdge,yEdge,normal(:,1),normal(:,2))
hold off

    
    

    


