function [root, residual, errStatus] = SecondOrderNewtonRaphsonRoot(fun, initGuess, relativeTol, scale, nMaxIter, failCriteria)
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

% Enables debugging lines, 0 is for disabled.
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
    nMaxIter = 1000;
end

% setup default fail criteria
if nargin <6
    failCriteria = @(x) false;
end

% Set up initial guess
x = initGuess;

% Compute the function and its 1st and 2nd derivatives
[f, fx, fxx] = Compute_f_fx_fxx(fun, x, dx);
fOld = 2*f;

% 2nd Order Zero Crossing condition
del = fx^2 - 2 * f * fxx;

% Compute step size required to reach 0
if fxx ~= 0
    if del >= 0
        xStep = -fx/fxx + sign( -f*fx*fxx) * sqrt(del)/fxx;
    else
        xStep = -f/fx * ( 1 + 0.5*f/fx^2*fxx );
    end
else
    
    xStep = -f/fx;
end

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
    %     error('Derivative is zero at this point. This is a local minimum.');
    % 2;
end

% Number of iterations
nIter = 0;

if DEBUG > 0
    % History
    xHist = zeros(1, nMaxIter);
    yHist = zeros(1, nMaxIter);
    xHist(nIter+1) = x;
    yHist(nIter+1) = f;
end

while abs(f/fOld-1) > relativeTol ...
        && errStatus == 0 ...
        && nIter < nMaxIter ...
        && failCriteria(x) == false
    
    % TODO if is evaluated at the end of previous iter, and start of this
    % iter. Can be optimized by inputting the value along with the function.
    % Compute the function and its 1st and 2nd derivatives
    [f, fx, fxx] = Compute_f_fx_fxx(fun, x, dx);
    
    % 2nd Order Zero Crossing condition
    del = fx^2 - 2 * f * fxx;
    
    % Compute step size required to reach 0
    if del >= 0
        % Evaluate both possible roots and choose the closest one
        delSqrt = sqrt(del);
        dx1 = -fx - delSqrt;
        dx2 = -fx + delSqrt;
        if fxx ~= 0
            if abs(dx1) < abs(dx2)
                xStep = dx1/fxx;
            else
                xStep = dx2/fxx;
            end
        else
            xStep = -f/fx;
        end
    else
        xStep = -f/fx * ( 1 + 0.5*f/fx^2*fxx );
    end
    
    if DEBUG > 1
        delXList = ( -100:100 ) * xStep/10;
        xList = x + delXList;
        
        %         ax^2 + bx + c
        a = fxx/2;
        b = fx;
        c = f;
        
        
        figure(2001);
        clf
        plot(xList, a*delXList.^2 + b*delXList + c, 'DisplayName','2nd Order Fit' );
        grid on
        title('2nd Order Fit')
        hold on
        plot( x*ones(2,1), ylim, 'r--','DisplayName', 'StartPoint' )
        plot( (x+xStep)*ones(2,1), ylim, 'g--', 'DisplayName','Target Step'  );
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
        %         errStatus = 4;
    end
    
    % Take underrelaxed step
    x = x + xStep * scale;
    
    % Store previous residual
    fOld = f;
    
    % Compute the value of fun
    f = fun(x);
    
    % Compute error status boundary exceeded
    if x < xMin || x > xMax
        errStatus = 2;
    end
    
    % Update Iteration Count
    nIter = nIter + 1;
    
    if DEBUG > 0
        % Record values for debugging
        xHist(nIter+1) = x;
        yHist(nIter+1) = f;
    end
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


function [f, fx, fxx] = Compute_f_fx_fxx(fun, x, dx)
% takes up to nth derivative of the function

% Center point
f = fun(x);

% Previous point
fm1 = fun(x-dx);

% Next Point
fp1 = fun(x+dx);

% Compute first derivative
fx = (fp1 - fm1)/(2*dx);

% Compute second derivative
fxx = (fp1+fm1-2*f)/dx^2;
end