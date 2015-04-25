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

% Get x and y points by finding nonzero elements
[xPoints, yPoints] = find(expandedDomain);

% Form operator matrix for equation
% Start with fours on diagonals
operatorMatrix = 4 * eye(numel(xPoints));

% Check surrounding points and use -1 where necessary
for i=1:numel(xPoints)
    currentPoint = [xPoints(i), yPoints(i)];
    
    % Set up coordinates of surrounding points
    pointAbove = [currentPoint(1) + 1, currentPoint(2)];
    pointBelow = [currentPoint(1) - 1, currentPoint(2)];
    pointLeft = [currentPoint(1), currentPoint(2) - 1];
    pointRight = [currentPoint(1), currentPoint(2) + 1];
    
    % Vector of surrounding points for easy iteration
    surroundingPoints = [pointAbove; pointBelow; pointLeft; pointRight];
    
    % Get single indices of surrounding points
    for j=1:numel(surroundingPoints) / 2
        idx = findPoint(xPoints, yPoints, [surroundingPoints(j, 1); surroundingPoints(j,2)]);
        operatorMatrix(i, idx) = -1;
    end
end

% Create RHS vector
rhsVector = rhs(xPoints, yPoints);

% Solve System
solutionMatrix = operatorMatrix \ rhsVector;

% Coordinates for Return
[X, Y] = meshgrid(stepSize:stepSize:size(domainMatrix, 2), ...
    stepSize:stepSize:size(domainMatrix, 1));

Z = expandedDomain;
for i=1:numel(solutionMatrix)
    Z(find(Z == 1, 1)) = solutionMatrix(i);
end
end

