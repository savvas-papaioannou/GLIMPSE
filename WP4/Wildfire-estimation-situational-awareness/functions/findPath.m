function path = findPath(map,covermap,dronePosition,resolution)
        
    covermap = covermap - (map & covermap);

    G = occupancyMatrixToGraph(map, resolution);

    
    [rows, cols] = find(covermap);
    nodesToVisit = sub2ind(size(covermap),cols,rows);
    
    [~,startIndex] = min(sqrt((cols*resolution - dronePosition(1)).^2 + (rows*resolution - dronePosition(2)).^2));

    numNodes = numel(nodesToVisit);
    shortestPaths = cell(numNodes); % Store paths here
    shortestPathLengths = inf(numNodes); % Store lengths here
    
    for i = 1:numNodes
        for j = 1:numNodes
            if i ~= j
                [path, d] = shortestpath(G, nodesToVisit(i), nodesToVisit(j));
                shortestPaths{i, j} = path;
                shortestPathLengths(i, j) = d;
            end
        end
    end
    
    
    % Greedy TSP solution - Start from the first node, then always go to the nearest unvisited node
    currentNode = startIndex; % Start from the calculated startIndex
    visited = false(1, numNodes);
    visited(currentNode) = true;
    orderOfVisit = currentNode;
    
    while sum(visited) < numNodes
        distancesToUnvisited = shortestPathLengths(currentNode,:);
        distancesToUnvisited(visited) = inf; % Ignore visited nodes
        [minDist, nextNode] = min(distancesToUnvisited);
        visited(nextNode) = true;
        orderOfVisit = [orderOfVisit nextNode];
        currentNode = nextNode;
    end
    
    % Retrieve the actual paths and total distance
    totalDistance = 0;
    finalPath = [];
    
    for i = 1:length(orderOfVisit)-1
        finalPath = [finalPath, shortestPaths{orderOfVisit(i), orderOfVisit(i+1)}(1:end-1)];
        totalDistance = totalDistance + shortestPathLengths(orderOfVisit(i), orderOfVisit(i+1));
    end
    finalPath = [finalPath, nodesToVisit(orderOfVisit(end))]; % Add the last node
    [X,Y] = ind2sub(size(covermap),finalPath);
    path = [X;Y]';

    %{
    figure(9);
    hold on
    [rows, cols] = find(covermap);
    % Create a scatter plot of these indices
    scatter(cols, rows, 'filled'); % Use 'filled' for filled markers
    [rows, cols] = find(map == 1);
    % Create a scatter plot of these indices
    scatter(cols, rows, 'filled'); % Use 'filled' for filled markers
    
    [rows, cols] = ind2sub(size(covermap),finalPath);
    plot(rows, cols, '-o'); % Use 'filled' for filled markers

    hold off
    xlim([1, size(covermap,1)]); % Adjust x-axis limits to match matrix size
    ylim([1, size(covermap,1)]); % Adjust y-axis limits to match matrix size
    % Adding labels and title for clarity
    xlabel('X-axis');
    ylabel('Y-axis');
    zlabel('Value');

    %% Plot
    % Generate XData and YData for plotting
    [rows, cols] = size(map);
    [Y, X] = meshgrid(1:rows,1:cols);
    XData = X(:)';
    YData = Y(:)'; % Negate YData to flip the plot to match the matrix orientation

    % Plot the graph
    figure(10);
    p = plot(G, 'XData', XData, 'YData', YData); % Plot the graph with specified node positions
    axis equal;
    title('Path Planning with Calculated Path Highlighted');
    xlabel('Column Index');
    ylabel('Row Index'); % Label adjusted for flipped Y-axis
    
    % Highlight the calculated path if it exists
    % Assuming 'path' is calculated from shortestpath(G, startNode, endNode)
    if ~isempty(path)
        highlight(p, finalPath, 'EdgeColor', 'r', 'LineWidth', 2);
    end
    %}
end


function G = occupancyMatrixToGraph(matrix, resolution)
    [rows, cols] = size(matrix);
    G = digraph;
    G = addnode(G, rows*cols);

    for  r = 1:rows
        for c = 1:cols
            if matrix(r, c) == 0 % Only proceed if current cell is not an obstacle
                currentNode = (r-1) * cols + c;
                for dr = -1:1
                    for dc = -1:1
                        if dr == 0 && dc == 0
                            continue; % Skip the current cell
                        end
                        newR = r + dr;
                        newC = c + dc;
                        if newR > 0 && newR <= rows && newC > 0 && newC <= cols && matrix(newR, newC) == 0
                            newNode = (newR-1) * cols + newC;
                            distance = resolution;
                            if dr ~= 0 && dc ~= 0
                                distance = sqrt(2) * resolution;
                            end
                            G = addedge(G, currentNode, newNode, distance);
                        end
                    end
                end
            end
        end
    end
end