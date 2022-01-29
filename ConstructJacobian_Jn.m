function J = ConstructJacobian_Jn(c,d,rho)
% Note that c(1) is c2
%           d(1) is d2

% Compute number of cols
nCols = size(rho, 2);

% Compute number of rows
nRows = 3 * (nCols- 2);

% Initialize Jacobian
J = zeros(nRows, nCols, 'like', rho);

% Fill for each 3*row
for i = 1 : nCols-2
    J( 3 * ( i-1 )+1 : 3 * i, i + (0:2) ) = [ c(i) * rho(:, i), - rho(:, i+1), d(i) * rho(:, i+2)];
end
end