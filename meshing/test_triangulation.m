clear;clc;close all

n = 200;
theta = linspace(0,2*pi,n+1);
theta(end) = [];

radius = 0.1+rand(size(theta)).^2;

x = radius.*cos(theta);
y = radius.*sin(theta);

% x = [...
%      0.060647
%        0.1195
%      0.027518
%    -0.0063566
%     -0.013888
%    -0.0045278
%      -0.52104
%     -0.064312
%       0.11168
%      0.089945];
%  
%  y = [...
%             0
%       0.08682
%      0.084691
%      0.019564
%       0.01009
%    5.5449e-19
%      -0.37856
%      -0.19793
%      -0.34372
%     -0.065349];


for i = 1:n
    text(x(i),y(i),num2str(i))
end

tic
i_tri = poly2tri(x,y);
toc
% hold on
% plot(x(i_tri),y(i_tri))

figure(898);clf
subplot(1,2,1)
patch('Faces',i_tri,'Vertices',[x;y]','facecolor','r','edgecolor','k')
hold on
plot(x,y)

%%
subplot(1,2,2)
% x = [0 1 0 1 0.5];    %# Unordered x coordinates of vertices
% y = [0 1 1 0 0.5];    %# Corresponding y coordinates of vertices
% edgeLines = [1 3;...  %# Point 1 connects to point 3
%              1 4;...  %# Point 1 connects to point 4
%              2 3;...  %# Point 2 connects to point 3
%              2 5;...  %# Point 2 connects to point 5
%              5 4];    %# Point 5 connects to point 4

edgeLines = [1:n;[2:n,1]]';

tic
dt = DelaunayTri(x(:),y(:),edgeLines);  %# Create a constrained triangulation
toc
isInside = inOutStatus(dt);  %# Find the indices of inside triangles
faces = dt(isInside,:);      %# Get the face indices of the inside triangles
vertices = [x(:) y(:)];      %# Vertex data for polygon
hPolygon = patch('Faces',faces,...
                 'Vertices',vertices,...
                 'FaceColor','r');  %# Plot the triangular faces in red

set(hPolygon,'EdgeColor','k');  %# Turn off the edge coloring
xEdge = x(edgeLines).';           %'# Create x coordinates for the edge
yEdge = y(edgeLines).';           %'# Create y coordinates for the edge
hold on;                           %# Add to the existing plot
line(xEdge,yEdge,'Color','k');     %# Plot the edge in black