function [X, Y, Z] = centralDifferencePoisson(stepSize, domainMatrix, rhs)
% Solve poisson equation with domain given by unit squares in
% matrix. The differential equation is then represented by a
% matrix system of equations to solve.

% Number of steps
steps = 1 / stepSize;

% Maximum x and y of domain for matrix dimensions
domainXMax = size(domainMatrix, 2);
domainYMax = size(domainMatrix, 1);

% Domain expanded by mesh
expandedDomain = zeros(domainYMax * steps, domainXMax * steps);

% Fill in expanded domain
inRow = 0;
for row = 1:domainYMax                  % Iterate thorough domain rows
    ones = 0;                           % Number of ones in this row
    start = -1;                         % Column of first one
    for col = 1:domainXMax
        if domainMatrix(row, col) == 1  % Check if element is 1
        ones = ones + 1;                % Increment ones
            if start == -1              % Set column if not yet set
            start = col;
            end
        end
    end
    if ones > 0                         % Add to expandedDomain if needed
        for step = 0:steps - 1          % Each step contributes to a row
            inRow = inRow + 1;          % Increment input row
            for inCol = 0:steps * ones - 1  % Input column
                expandedDomain(inRow, start + inCol) = 1;
            end
        end
    end
end

% Get mesh points
% A mesh point occurs at the intersection of four boxes,
% which are represented as 1s in expandedDomain and domainMatrix.
% The coordinates of these are given by the top right square, which
% is given by the row of the square plus the corresponding row of the
% unit mesh and the column of the square plus the corresponding column
% of the unit mesh

% Rows and columns of expandedDomain for easy access
expandedDomainRows = size(expandedDomain, 1);
expandedDomainCols = size(expandedDomain, 2);

% Index of mesh point. Start at 0 since we have 0 mesh points at the
% beginning
meshPointIndex = 0;
% Empty list for mesh points
meshPoints = [];
    
% Iterate through rows backwards due to the way the matrices are indexed
for row = expandedDomainRows:-1:2
    % Iterate through columns
    for col = 1:expandedDomainCols - 1
        % Check if top right square is 1. If so, this mesh point is given
        % by the coordinates of that square scaled by the step size
        if expandedDomain(row, col) == 1 && ...
                expandedDomain(row - 1, col + 1) == 1
            meshPointIndex = meshPointIndex + 1;
            meshPoints(meshPointIndex,1) = stepSize * col;
            meshPoints(meshPointIndex,2) = ...
                stepSize * (expandedDomainRows - row + 1);
        end
    end
end

% Form operator matrix for equation
% Start with fours on diagonals
operatorMatrix = 4 * eye(meshPointIndex);

% Assign negative 1s to opperator matrix
for row = 1:meshPointIndex
    for col = 1:meshPointIndex
        % Check if there are any adjacent points within range
        if (meshPoints(row, 1) == meshPoints(col, 1) + stepSize ...
                && meshPoints(row, 2) == meshPoints(col, 2)) || ...
                (meshPoints(row, 1) == meshPoints(col, 1) - stepSize ...
                && meshPoints(row, 2) == meshPoints(col, 2)) || ...
                (meshPoints(row, 1) == meshPoints(col, 1) ...
                && meshPoints(row, 2) == meshPoints(col, 2) + stepSize) || ...
                (meshPoints(row, 1) == meshPoints(col, 1) ...
                && meshPoints(row, 2) == meshPoints(col, 2) - stepSize)
            operatorMatrix(row, col) = -1;
        end
    end
end

% Create RHS vector
rhsVector = [];
for index = 1:meshPointIndex
    rhsVector(index) = stepSize^2 * ...
        rhs(meshPoints(index, 1) ,meshPoints(index, 2));
end

% Solve System
solutionMatrix = operatorMatrix \ transpose(rhsVector);

% Coordinates for Return
[X, Y] = meshgrid(0:stepSize:size(domainMatrix, 2), ...
    0:stepSize:size(domainMatrix, 1));
Z = griddata(meshPoints(:,1), meshPoints(:,2), solutionMatrix, X, Y);
% TODO: Spy operator matrix
% TODO: Fix domain step (maybe griddata?)
end

