function [root, residual, errStatus] = NewtonRaphsonRoot(fun,dfun, initGuess, tol, scale, nMaxIter)
%NewtonRaphsonRoot Find the root of a function given its derivative
%function
%   Inputs:
%       fun             : function handle
%       dfun            : function handle for the derivative of the function
%       initGuess       : initial guess
%       tol(optional)   : tolerance for the final result (default 1E-5)
%       scale(optional) : relaxation scale (default 0.5)
%
%   Outputs:
%       root    : closest root to the initial guess
%       residual: absolute value of residual of the root (if negative, then error exists)

DEBUG = true;

% setup tolerance if not given
if nargin < 4
    tol = 1E-5;
end

% setup under relaxation if not given
if nargin < 5
    scale = 0.5;
end

% setup maximum number of iterations
if nargin < 6
    nMaxIter = 100;
end

% Number of maximum consecutive oscillations of the derivative
nMaxOsc = 5;

% Set up initial guess
x = initGuess;

% Compute the value of fun
y = fun(x);

% Compute the derivative of the function at the point
dydxOld = dfun(x);

% Compute step size required to reach 0
xStep = -y/dydxOld;
    
% Error Status
errStatus = 0;

% Determine if x is lower or upper bound
if xStep > 0
    xMin = x;
    xMax = Inf;
elseif xStep < 0
    xMin = -Inf;
    xMax = x;
else
    error('Derivative is zero at this point. This is a local minimum.');
end

% Number of iterations
nIter = 0;

% Number of consecutive times derivative changed signs;
nOsc = 0;

% History
xHist = zeros(1, nMaxIter);
yHist = zeros(1, nMaxIter);
xHist(nIter+1) = x;
yHist(nIter+1) = y;

while abs(y) > tol && errStatus == 0 && nIter < nMaxIter

    % Compute the derivative of the function at the point
    dydx = dfun(x);

    % Compute step size required to reach 0
    xStep = -y/dydx;
    
    % Determine if x is lower or upper bound
    if xStep > 0
        xMin = x;
    elseif xStep < 0
        xMax = x;
    else
        error('Step is zero.') % This should not happen as abs(y) > tol
    end

    % Take underrelaxed step
    x = x + xStep * scale;

    % Compute the value of fun
    y = fun(x);
    
    % Compute consecutive times derivative changed signs
    if dydxOld * dydx < 0
        nOsc = nOsc + 1;
    else
        nOsc = 0;
    end
    
    
    % Compute error status 5 oscillations
    if nOsc > nMaxOsc
        errStatus = 1;
    end
    
    % Compute error status boundary exceeded
    if x < xMin || x > xMax
        errStatus = 2;
    end
    
    % Update previous derivative
    dydxOld = dydx;
    
    % Update Iteration Count
    nIter = nIter + 1;

    % Record values for debugging
    xHist(nIter+1) = x;
    yHist(nIter+1) = y;
end

% Compute error status maximum iterations exceeded
if nIter >= nMaxIter
    errStatus = 3;
end
    
% Set outputs
root = x;
residual = abs(y);

% Set residual to negative for error
if errStatus ~= 0
    residual = -residual;
end

if DEBUG
    figure(1001);
    clf
    plot(1:nIter+1, xHist(1:nIter+1))
    title('Root Guess with iterations')
    xlabel('No Iter')
    grid on
    
    figure(1002);
    clf
    plot(1:nIter+1, yHist(1:nIter+1))
    title('Residual with iterations')
    xlabel('No Iter')
    grid on
    
    figure(1003);
    clf
    plot(xHist(1:nIter+1), yHist(1:nIter+1))
    title('Function')
    grid on
    
end

end

