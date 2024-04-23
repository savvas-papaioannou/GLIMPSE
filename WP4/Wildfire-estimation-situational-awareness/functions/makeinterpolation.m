
function [interpolatedVector] = makeinterpolation(largerVector,smallerVector)


% Interpolate the smaller vector to match the size of the larger one
interpSize = size(largerVector, 1);

% Initialize arrays to hold the interpolated points
interpolatedSmallerVector = NaN(interpSize, size(smallerVector, 2));

% Original positions in the interpolated vector
originalPositions = linspace(1, interpSize, size(smallerVector, 1));

% Loop through each dimension separately for interpolation
for dim = 1:size(smallerVector, 2)
    % Pre-fill with NaNs
    interpolatedColumn = NaN(interpSize, 1);

    % Place original values at corresponding positions
    interpolatedColumn(round(originalPositions)) = smallerVector(:, dim);

    % Find NaN indices for interpolation
    nanIndices = isnan(interpolatedColumn);

    % Generate x values for existing points
    xExisting = find(~nanIndices);

    % Interpolate only NaN values
    interpolatedColumn(nanIndices) = interp1(xExisting, interpolatedColumn(~nanIndices), find(nanIndices), 'linear');

    % Assign interpolated column to the output matrix
    interpolatedSmallerVector(:, dim) = interpolatedColumn;
end

% Extract the first element of the larger vector
firstElementLargerVector = largerVector(1, :);

% Calculate Euclidean distances
distances = sqrt(sum((interpolatedSmallerVector - firstElementLargerVector).^2, 2));

% Find the index of the closest point
[~, closestIndex] = min(distances);

% Reorder the interpolatedSmallerVector to make the closest point the first element
% and preserve the cyclical nature
if closestIndex == 1
    % No need to reorder if the closest point is already the first element
    reorderedVector = interpolatedSmallerVector;
else
    % Otherwise, reorder taking into account the cyclical nature
    reorderedVector = [interpolatedSmallerVector(closestIndex:end-1, :); interpolatedSmallerVector(1:closestIndex-1, :); interpolatedSmallerVector(closestIndex, :)];
end

interpolatedVector = reorderedVector;

%interpolatedVector = interpolatedSmallerVector;