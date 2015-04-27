function [ X, Y, Z ] = explicit_h1_solver( rhs )
%EXPLICITH1_SOLVER Explicitly solves the poisson equation on domain
%   Explicitly solves the poisson equation on the given domain with step
%   size of 1. rhs Is a function pointer with two arguments

% Differential Operator Matrix
D = 4 * eye(7);
D(1,2) = -1;
D(1,4) = -1;
D(2,1) = -1;
D(2,5) = -1;
D(2,3) = -1;
D(3,2) = -1;
D(4,1) = -1;
D(4,5) = -1;
D(6,4) = -1;
D(7,6) = -1;

% X and Y Coordinates
X = [1; 2; 3; 4; 1; 2; 1];
Y = [1; 1; 1; 1; 2; 2; 3];

% Evaluate rhs at points
rhsVector = rhs(X,Y);

% Solve for solution vector
solutionVector = D \ rhsVector;

[X,Y] = meshgrid(0:5, 0:4);

Z = [ 0, 0, 0, 0, 0, 0; 
    0, solutionVector(1), solutionVector(2), solutionVector(3), solutionVector(4), 0;
    0, solutionVector(5), solutionVector(6), 0, 0, 0;
    0, solutionVector(7), 0, 0, 0, 0;
    0, 0, 0, 0, 0, 0];

% Get rid of x and y coordinates not inside domain or boundary
zeroIndices = find(~conv2(Z, [1 1 1; 1 0 1; 1 1 1], 'same'));
X(zeroIndices) = NaN;
Y(zeroIndices) = NaN;

end