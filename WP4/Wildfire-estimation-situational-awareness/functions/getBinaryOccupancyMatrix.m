function occupancyGrid = getBinaryOccupancyMatrix(points,gridSizeX,gridSizeY,gridResolution)
    %% Calculates Binary Occupancy Matrix
    % Inputs:
    % points - fire perimeter points Nx2 matrix where the first column is x and the second is y.
    % grid parameters (size and resolution)
    % gridSizeX - x coordinate size in meters
    % gridSizeY - y coordinate size in meters
    % gridResolution - size of each grid cell (assuming meters)
    
    % Outputs:
    % occupancyGrid - binary occupancy grid, each cell inside the perimeter is 1
    % and outsite 0
    
    % Create an empty occupancy grid with the correct dimensions
    occupancyGrid = zeros(ceil(gridSizeX / gridResolution), ceil(gridSizeY / gridResolution));
    
    % Scale points according to the grid resolution
    scaledPointsX = ceil(points(:, 1) / gridResolution);
    scaledPointsY = ceil(points(:, 2) / gridResolution);
    
    % Define the grid of points to test
    [xGrid, yGrid] = meshgrid(1:size(occupancyGrid, 2), 1:size(occupancyGrid, 1));
    
    % Use inpolygon to find points inside the scaled boundary
    inside = inpolygon(xGrid, yGrid, scaledPointsX, scaledPointsY);
    
    % Fill the occupancy grid
    occupancyGrid(inside) = 1;
end