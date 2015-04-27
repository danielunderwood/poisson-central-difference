function [X, Y, Z] = centralDifferencePoisson(stepSize, domainMatrix, rhs)
% Solve poisson equation with domain given by unit squares in
% matrix. The differential equation is then represented by a
% matrix system of equations to solve.

% Number of steps
steps = round(1 / stepSize);

% Expand domain using the kronecker tensor product
expandedDomain = kron(domainMatrix, ones(steps));

% Fix x axis missing expanded values and boundaries
leftSide = kron(expandedDomain(:,1), ones(1, steps-1));
leftBoundary = zeros(size(leftSide,1), 1);
rightBoundary = leftBoundary;
expandedDomain = [leftBoundary, leftSide, expandedDomain, rightBoundary];

% Fix y axis missing expanded values and boundaries
bottom = kron(expandedDomain(1,:), ones(steps-1, 1));
topBoundary = zeros(1, size(bottom,2));
bottomBoundary = topBoundary;
expandedDomain = [bottomBoundary; bottom; expandedDomain; topBoundary];

% Get x and y points by finding nonzero elements
[yPoints, xPoints] = find(expandedDomain);

% Scale x and y points based on step size
xPoints = stepSize * xPoints;
yPoints = stepSize * yPoints;

% Form operator matrix for equation
% Start with fours on diagonals
operatorMatrix = 4 * eye(numel(xPoints));

% Check surrounding points and use -1 where necessary
for i=1:numel(xPoints)
    currentPoint = [xPoints(i), yPoints(i)];
    
    % Set up coordinates of surrounding points
    pointAbove = [currentPoint(1), currentPoint(2) + stepSize];
    pointBelow = [currentPoint(1), currentPoint(2) - stepSize];
    pointLeft = [currentPoint(1) - stepSize, currentPoint(2)];
    pointRight = [currentPoint(1) + stepSize, currentPoint(2)];
    
    % Vector of surrounding points for easy iteration
    surroundingPoints = [pointAbove; pointBelow; pointLeft; pointRight];
    
    % Get single indices of surrounding points
    for j=1:numel(surroundingPoints) / 2
        idx = findPoint(xPoints, yPoints, [surroundingPoints(j, 1); surroundingPoints(j,2)]);
        operatorMatrix(i, idx) = -1;
    end
end

% Create RHS vector
rhsVector = stepSize.^2 * rhs(xPoints, yPoints);

% Solve System
solutionMatrix = operatorMatrix \ rhsVector;

% Coordinates for Return
[X, Y] = meshgrid(0:stepSize:size(domainMatrix, 2)+1, ...
    0:stepSize:size(domainMatrix, 1)+1);

Z = expandedDomain;
for i=1:numel(solutionMatrix)
    Z(find(Z == 1, 1)) = solutionMatrix(i);
end

% Get rid of x and y coordinates not inside domain or boundary
zeroIndices = find(~conv2(Z, [1 1 1; 1 0 1; 1 1 1], 'same'));
X(zeroIndices) = NaN;
Y(zeroIndices) = NaN;
    

end

