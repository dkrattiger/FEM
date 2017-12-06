function i_triangle = poly2tri(x,y)

% use the ear-clipping method of triangulating polygons (this is quite
% simple and works on convex polygons)

% main idea:
% Every simple polygon (no self intersections or holes) contains at least
% one "ear". An ear is defined as a triangle with two edges along the
% polygon edge, and the last edge entirely contained within the polygon.
%
% algorithm:
% (1) select a convex vertex and check if the corresponding ear contains any
% internal points.
% 
%   (a) if so, it is not an ear, select a new convex vertex and repeat
%   (b) if not, it is an ear. Remove the ear and start over until only
%   three vertices remain


n = length(x);

% first check if polygon is CCW
CCW_poly = isccw(x,y);

%% create a list of convex vertices
convlist = zeros(1,n);
for i = 1:n
    index3 = rem([i-1,i,i+1]+n-1,n)+1;
    convlist(i) = (isccw(x(index3),y(index3)) == CCW_poly);
end

i_rem = 1:n;
n_rem = n;
i_triangle = zeros(n-2,3);



% populate indices for n-2 triangles
for i = 1:(n-2)  

    % loop over vertices in remaining polygon to find an "ear"
    for j = 1:n_rem 
        
        % if the jth node is convex...
        if convlist(i_rem(j))
            
            % form index for triplet of points (using points in remaining
            % polygon)
            index3 = i_rem(rem([j-1,j,j+1]+n_rem-1,n_rem)+1);
            others = i_rem; 
            others(rem([j-1,j,j+1]+n_rem-1,n_rem)+1) = [];
            
            % if no points are in triangle...
            if ~any(intriangle(x(index3),y(index3),x(others),y(others)))
                
                
                % add the corresponding triangle to triangle index list
                i_triangle(i,:) = index3;
                
                % remove ear for jth remaining vertex from remainder
                % polygon
                i_rem(j) = [];
                n_rem = n_rem - 1;
                
                % update convexity list (set removed node to false)
                convlist(index3(2)) = false;
                
                % update convexity of previous vertex
                index3 = i_rem(rem([j-2,j-1,j]+n_rem-1,n_rem)+1);
                convlist(index3(2)) = (isccw(x(index3),y(index3)) == CCW_poly);
                
                % update convexity of previous vertex
                index3 = i_rem(rem([j-1,j,j+1]+n_rem-1,n_rem)+1);
                convlist(index3(2)) = (isccw(x(index3),y(index3)) == CCW_poly);
                break
            end
        end
    end
end
   
function out = isccw(x,y)
    
    % length of point list
    n = length(x);

    % Sum of (x2-x1)*(y2+y1) will be positive for CW
    testval = sum((x([2:n,1])-x([1:n])).*(y([2:n,1])+y([1:n])));
    if  testval > 0
        out = false;
    else
        out = true;
    end

function result = intriangle(xtri,ytri,xp,yp)

    
    % make sure that if a list of points is given to check, that the list
    % is a vertical vector
    xp  = xp(:);
    yp  = yp(:);
    
    % initialize output vector
    result = true(size(xp));
    for i = 1:3
        
        index = rem([i-1,i,i+1]+3-1,3)+1;
        
        % create vectors from current vertex to other triangle vertices,
        % va and vb
        va = [xtri(index(3))-xtri(index(2));ytri(index(3))-ytri(index(2))];
        vb = [xtri(index(1))-xtri(index(2));ytri(index(1))-ytri(index(2))];
        
        % create vectors from current vertex to all test points
        vp = [xp.'-xtri(index(2));yp.'-ytri(index(2))];
        
        % compute the cross product between va and vb
        cross1 = va(1)*vb(2) - va(2)*vb(1);
        
        % compute the cross products between va and all test point vectors
        cross2 = va(1)*vp(2,:) - va(2)*vp(1,:);
        
        % if cross-products both have same sign, then the point is on the
        % correct side of va
        
        test = cross2/cross1 > 0;
        
        % update test with in/out value
        result = result & test;
    end