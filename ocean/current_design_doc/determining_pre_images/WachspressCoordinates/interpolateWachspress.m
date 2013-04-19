clear all

%interpolate with Wachspress coordinates
%given a hexagon, compare reconstruction of test function to exact soln

% when test function is linear in x and y
% reconstruction should be at the level of round off error

% the routine plot the error at the of the computation

nX=50;
nY=50;
reconstruct=NaN(nX,nY);
eval=NaN(nX,nY);
diff=NaN(nX,nY);

nTri=3 % number of vertices on a triangle
nDim=2 % number of spatial dimensions

%build polygon (here it is a perfect hexagon)
nVertices=6
b = sqrt(3)/2;
x=transpose(([-1/2 1/2 1 1/2 -1/2 -1]))
y=transpose(([-b -b 0 b b 0]))

%the polygon is broken into triangles
%the number of triangles is nVertices
tri = zeros(nVertices,nTri,nDim);

%find the B's of the Wachspress coordinate
for i=1:nVertices
    iP1=i+1;
    iM1=i-1;
    if i==nVertices
        iP1=1;
    end
    if i==1
        iM1=nVertices;
    end
    
    tri(1,1,1)=x(i);
    tri(1,1,2)=y(i);
    tri(1,2,1)=x(iP1);
    tri(1,2,2)=y(iP1);
    tri(1,3,1)=x(iM1);
    tri(1,3,2)=y(iM1);
    
    work(:,:) = tri(1,:,:);
    [ B(i,1) ] = areaBasedonVertexCoords( work );
    
end

%evaluate (linear) function at each vertex point
for i=1:nVertices
    [ q(i,1) ] = testFunction( x(i,1), y(i,1) );  
end

%scan across the polygon
%if the coordinate is inside polygon, then reconstruct function
%compare the reconstruction to the evaluation of the test function
for jP=1:nY
    for iP=1:nX
        
        yP=-b + b*(jP-1)/((nY-1)/2);
        xP=-1.0+(iP-1)/((nX-1)/2);
        
        rFrac = abs(yP)/b;
        xEdge = 1.0 - rFrac/2.0;
        
        if abs(xP)<xEdge
        
        %find triangle area, tri(nTri,nVert,x)
        for i=1:nVertices
          iP1=i+1;
          if i==nVertices
              iP1=1;
          end
          tri(i,1,1)=xP;
          tri(i,1,2)=yP;
          tri(i,2,1)=x(i,1);
          tri(i,2,2)=y(i,1);
          tri(i,3,1)=x(iP1,1);
          tri(i,3,2)=y(iP1,1);
    
          % find area of triangles (the A's of Wachspress coordinates)
          work(:,:) = tri(i,:,:);
          [ A(i,1) ] = areaBasedonVertexCoords( work );
    
        end

        %find the unnormalized wachspress weights
        for i=1:nVertices
    
         iM1=i-1;
         if i==1
           iM1=nVertices;
         end
   
         AMask = A;
         AMask(i)=1.0;
         AMask(iM1)=1.0;
    
         w(i,1)=1.0;
         for j=1:nVertices
             w(i,1) = w(i,1)*AMask(j,1);
         end
    
         w(i,1) = B(i,1)*w(i,1);
        end

        total = sum(w(:,1));
        gamma(:,1) = w(:,1) / total;
        
        reconstruct(iP,jP)=0.0;
        for i=1:nVertices
            reconstruct(iP,jP) = reconstruct(iP,jP) + gamma(i,1)*q(i,1);
        end
        
        [ eval(iP,jP) ] = testFunction( xP, yP );
        
        end
        
    end
end

diff = reconstruct - eval;
contourf(diff)
colorbar