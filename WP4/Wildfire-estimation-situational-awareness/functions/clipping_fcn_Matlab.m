function pNoHoles = clipping_fcn_Matlab(q_in)
    % Create the initial polyshape without enforcing boundary orientation
    pout = simplify(polyshape(q_in(:,1), q_in(:,2), "Simplify",false));
    prmholes = rmholes(pout);
    
    % Sort the regions by area in descending order
    pSorted = sortregions(prmholes,'area','descend');
    R = regions(pSorted);
    
    % Extract the largest region's vertices
    largestRegionVertices = R(1).Vertices;
    
    % Check and adjust the orientation to counter-clockwise if necessary
    if isClockwise(largestRegionVertices)
        largestRegionVertices = flipud(largestRegionVertices);
    end
    
    % Create the polyshape with the possibly adjusted vertices
    pccw  = polyshape(largestRegionVertices(:,1), largestRegionVertices(:,2), 'Simplify', false);
    pNoHoles = [pccw.Vertices;pccw.Vertices(1,:)];
end

% Helper function to determine if the vertices are ordered clockwise
function cw = isClockwise(vertices)
    area = 0;
    n = size(vertices, 1);
    for i = 1:n-1
        area = area + (vertices(i,1) * vertices(i+1,2) - vertices(i+1,1) * vertices(i,2));
    end
    area = area + (vertices(n,1) * vertices(1,2) - vertices(1,1) * vertices(n,2)); % Close the polygon
    area = area / 2;
    cw = area < 0;
end