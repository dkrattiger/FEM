function iTri = patch3D2tri(varargin)

% This function performs a triangulation of a 3D patch. It uses constrained
% Delaunay triangulation to perform the triangulation. This only works on
% 2D polynomials so the points are first mapped to a plane
%
% iTri = patch3Dtri(X,Y,Z)
% iTri = patch3Dtri(C)      (C = [x,y,z])

if nargin == 3
    x = varargin{1};
    y = varargin{2};
    z = varargin{3};
    x = x(:); y = y(:); z = z(:);
    X = [x,y,z];
    n = length(x);
end

if nargin == 1
    X = varargin{1};
    [n,m] = shape(X);
    if n < m
        X = X';
        [n,m] = shape(X);
    end
end

% create average normal vector for patchface
N = sum(cross(X-X([n,1:n-1],:),X([2:n,1],:)-X),1)';
N = N/norm(N);

% create basis vectors orthogonal to normal vector
V1 = rand(3,1);
V1 = V1-N*(V1'*N);
V1 = V1/norm(V1);

V2 = cross(N,V1);
V2 = V2/norm(V2);

% find coordinates in V basis
% size(X)
% size(V1)
% size(V2)
% [V1,V2,N]*[z1,z2,z3,..] = [x1,x2,..]


Z = X/[V1,V2]';

% use constrained delaunay triangulation to triangulate the 2D polynomial 
% given by coordinates Z

% constraint lines (lines that are forced into triangulation)
edgelines = [1:n;[2:n,1]]';
DT = delaunayTriangulation(Z(:,1),Z(:,2),edgelines);

% pick out interior triangles
IO = isInterior(DT);

% pass connectivity list to return
iTri = DT.ConnectivityList(IO,:);