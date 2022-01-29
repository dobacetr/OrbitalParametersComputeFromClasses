function [root, residual, errStatus] = SecondOrderNewtonRaphsonRoot(fun, dfun, d2fun, initGuess, tol, scale, nMaxIter)
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

DEBUG = 2;

% setup tolerance if not given
if nargin < 5
    tol = 1E-5;
end

% setup under relaxation if not given
if nargin < 6
    scale = 0.5;
end

% setup maximum number of iterations
if nargin < 7
    nMaxIter = 100;
end

% Set up initial guess
x = initGuess;

% Compute the value of fun
f = fun(x);

% Compute the derivative of the function at the point
fx    = dfun(x);
fxx  = d2fun(x);

% Compute step size required to reach 0
xStep = -f/fx * ( 1 + 0.5*f/fx^2*fxx );
    
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


% History
xHist = zeros(1, nMaxIter);
yHist = zeros(1, nMaxIter);
xHist(nIter+1) = x;
yHist(nIter+1) = f;

while abs(f) > tol && errStatus == 0 && nIter < nMaxIter
    % Compute the value of fun
    f = fun(x);

    % Compute the derivative of the function at the point
    fx    = dfun(x);
    fxx  = d2fun(x);

    % Compute step size required to reach 0
    xStep = -f/fx * ( 1 + 0.5*f/fx^2*fxx );
    
    if DEBUG > 1
        xList = x + ( -100:100 ) * xStep/10;
        
        figure(2001);
        clf
        plot(xList, fxx/2*xList.^2 + fx*xList + f - (fxx/2*x^2+fx*x) );
        grid on
        title('2nd Order Fit')
        hold on
        plot( x*ones(2,1), ylim, 'r--','DisplayName', 'StartPoint' )
        plot( (x-xStep)*ones(2,1), ylim, 'g--', 'DisplayName', '1 Step'  );
        plot( (x+xStep)*ones(2,1), ylim, 'g--', 'HandleVisibility','off'  );
        legend Location best
    end
    
    % Determine if x is lower or upper bound
    if xStep > 0
        xMin = x;
    elseif xStep < 0
        xMax = x;
    else
        error('Step is zero.') % as abs(f) > tol, f/fx^2*fxx == -2 must be true.
        % TO DO: assign further meaning to this error
    end

    % Take underrelaxed step
    x = x + xStep * scale;

    % Compute the value of fun
    f = fun(x);
    
    % Compute error status boundary exceeded
    if x < xMin || x > xMax
        errStatus = 2;
    end
    
    % Update Iteration Count
    nIter = nIter + 1;

    % Record values for debugging
    xHist(nIter+1) = x;
    yHist(nIter+1) = f;
end

% Compute error status maximum iterations exceeded
if nIter >= nMaxIter
    errStatus = 3;
end
    
% Set outputs
root = x;
residual = abs(f);

% Set residual to negative for error
if errStatus ~= 0
    residual = -residual;
end

if DEBUG > 0
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

