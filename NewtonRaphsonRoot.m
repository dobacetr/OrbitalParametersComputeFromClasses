function [root, residual, errStatus] = NewtonRaphsonRoot(fun, initGuess, relativeTol, scale, nMaxIter, failCriteria)
%NewtonRaphsonRoot Find the root of a function given its derivative
%function
%   Inputs:
%       fun                     : function handle
%       initGuess               : initial guess
%       relativeTol(optional)   : tolerance of change ratio between
%                                 iterations abs(f/fold-1) (default 1E-5)
%       scale(optional)         : relaxation scale (default 0.5)
%       nMaxIter(optional)      : number of maximum iterations allowed (default 100)
%       failCriteria(optinal)   : a function handle that takes root as input and
%                                 returns logical. If true, execution will
%                                 stop. (default none)
%
%   Outputs:
%       root    : closest root to the initial guess
%       residual: absolute value of residual of the root (if negative, then error exists)
%       errStatus: 1 Failure Criteria is met
%                  2 root went out of scope
%                  3 Max Iter exceeded
%                  4 Function reached zero value


DEBUG = 0;

dx    = initGuess*1E-6;

% setup tolerance if not given
if nargin < 3
    relativeTol = 1E-5;
end

% setup under relaxation if not given
if nargin < 4
    scale = 0.5;
end

% setup maximum number of iterations
if nargin < 5
    nMaxIter = 100;
end

% setup default fail criteria
if nargin <6
    failCriteria = @(x) false;
end

% Set up initial guess
x = initGuess;

% Compute the function and its 1st derivative
f = fun(x);
fOld = 2*f;
fx = Compute_fx(fun, x, dx);

% Compute step size required to reach 0
xStep = -f/fx;
    
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

while abs(f/fOld-1) > relativeTol ...
        && errStatus == 0 ...
        && nIter < nMaxIter ...
        && failCriteria(x) == false
   
    % Compute the 1st derivative
    fx = Compute_fx(fun, x, dx);

    % Compute step size required to reach 0
    xStep = -f/fx;
    
    if DEBUG > 1
        xList = x + ( -100:100 ) * xStep/10;
        
        figure(2001);
        clf
        plot(xList, fx*xList + f - (fx*x) );
        grid on
        title('1st Order Fit')
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
%         error('Step is zero.') % as abs(f) > tol, f/fx^2*fxx == -2 must be true.
        % TO DO: assign further meaning to this error
        errStatus = 4; % Function became zero
    end

    % Take underrelaxed step
    x = x + xStep * scale;

    % Store old residual
    fOld = f;
    
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

% Computation has stopped due to exceeding fail criteria
if failCriteria(x) == true
    errStatus = 1;
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


function fx = Compute_fx(fun, x, dx)
% takes up to nth derivative of the function

% Previous point
fm1 = fun(x-dx);

% Next Point
fp1 = fun(x+dx);

% Compute first derivative
fx = (fp1 - fm1)/(2*dx);
end