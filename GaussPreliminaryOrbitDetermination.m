function [a, e, i, RAAN, w, TA, h, Eps, T] = GaussPreliminaryOrbitDetermination(time, rho, RSite, mu, flagIterationMode)

% Check the number of given observations
assert( length(rho) == 3 ...
        && length(time) == 3 ...
        , 'Only 3 Measurements can be implemented in the Gauss Method.' );

% Check the observation times, verify that they are in ascending order
assert(all(diff(time)>0), 'Measurement sequence must be ordered in ascending time.');
    
% 1) Compute Time differences
tau1 = time(1) - time(2);
tau3 = time(3) - time(2);
tau  = time(3) - time(1);

% 2) Compute Planes for observation pairs
P = zeros(3, 3);
P(:, 1) = cross(rho(:, 2), rho(:, 3));
P(:, 2) = cross(rho(:, 1), rho(:, 3));
P(:, 3) = cross(rho(:, 1), rho(:, 2));

% 3) Compute projection of 1st measurement on to the plane of 2nd and 3rd
D0 = dot(rho(:, 1), P(:, 1));
assert(abs(D0)>1E-12, 'Measurements can not be coplanar.');

% 4) Compute Dij
D = zeros(3,3);
for i = 1 : 3
    for j = 1 : 3
        D(i, j) = dot( RSite(:, i), P(:, j));
    end
end

% 5) Compute A and B
A = ( -D(1,2) * tau3/tau + D(2,2) + D(3,2)*tau1/tau )/D0;

B = ( D(1,2) * (tau3^2-tau^2)*tau3/tau + D(3,2)*(tau^2-tau1^2)*tau1/tau ) / (6*D0);

% 6) Compute E and Rsqr
E = dot(RSite(:, 2), rho(:, 2));

R2Sqr= dot( RSite(:, 2), RSite(:, 2) );

%  7) Compute a, b, c
a = - (A^2 + 2*A*E + R2Sqr);
b = -2*mu*B*(A+E);
c = -(mu*B)^2;

% 8) Compute r2
fun     = @(r2) r2.^8 + a*r2.^6 + b*r2.^3 + c; % = 0
initGuess = 6.4E6; % at about 700km altitude

% Try to find a real root
errStatus = 1;
r2 = 0;
% We dont check for real r2 because NewtonRaphson can not compute complex
% number
while ( ...
        errStatus ~= 0 ...      % Check for errors
        || r2 < 6378137 ...     % Check for positive altitude. This is already checked by the methods fail criteria
      ) ...
        && ...
      ( ...
        initGuess < 100E6 ...   % Check for initial guess attemps 
      )
    
%     [r2, residual, errStatus] = NewtonRaphsonRoot(fun, initGuess, 1E-8, 1.0, 100, @(x) x < 6378137);
    [r2, residual, errStatus] = SecondOrderNewtonRaphsonRoot(fun, initGuess, 1E-8, 1, 100, @(x) x < 6378137);
    initGuess = initGuess + 1E6;
end
% roots([1, 0, a, 0, 0, b, 0 , 0, c])
assert( initGuess < 100E6, 'Exceeded limit for initial guess to NewtonRapson.' );

% 9) compute measurement directions
rhoMag1 = ( ( 6 * ( D(3,1)*tau1/tau3 + D(2,1)*tau/tau3 )*r2^3 + mu*D(3,1)*( tau^2 - tau1^2 )*tau1/tau3 ) / ( 6 * r2^3 + mu * ( tau^2 - tau3^2 ) ) - D(1, 1) ) / D0;
rhoMag2 = A + mu*B/r2^3;
rhoMag3 = ( ( 6 * ( D(1,3)*tau3/tau1 - D(2,3)*tau/tau1 )*r2^3 + mu*D(1,3)*( tau^2 - tau3^2 )*tau3/tau1 ) / ( 6 * r2^3 + mu * ( tau^2 - tau1^2 ) ) - D(3, 3) ) / D0;

% 10) Compute r1, r2, r3
r1 = RSite(:, 1) + rhoMag1 * rho(:, 1);
r2 = RSite(:, 2) + rhoMag2 * rho(:, 2);
r3 = RSite(:, 3) + rhoMag3 * rho(:, 3);

% Iterative part
if flagIterationMode

    % 11) Lagrange coefficients Method
    r2Mag = sqrt(r2'*r2);
    f1 = 1 - 0.5 * mu / r2Mag^3 * tau1^2;
    f3 = 1 - 0.5 * mu / r2Mag^3 * tau3^2;
    g1 = tau1 - mu/( 6*r2Mag^3 ) * tau1^3;
    g3 = tau3 - mu/( 6*r2Mag^3 ) * tau3^3;

    v2 = ( -f3*r1 + f1*r3 ) / ( f1*g3 - f3*g1 );
        
    absRelativeError = Inf;
    
    nIter = 0;
    nMaxIter = 10000;
    tol = 1E-12;    % ratio
    
    oldRhos = [rhoMag1;rhoMag2;rhoMag3];
    
    while absRelativeError > tol ...
            && nIter < nMaxIter
        
        [a, e, i, RAAN, w, TA, h, Eps, T] = ComputeKepplerianElements(r2, v2, mu);
        if e >= 1
            error('eccentricity >= 1 is not accounted for.');
        end
    
        theta2 = TA;
        Me2 = ComputeMeanAnomalyFromTrueAnomaly(theta2, e);

        t1 = time(1);
        t2 = time(2);
        t3 = time(3);

        Me1 = Me2 + 2*pi*( t1 - t2 ) / T;
        Me3 = Me2 + 2*pi*( t3 - t2 ) / T;

        theta1 = ComputeTrueAnomalyFromMeanAnomaly(Me1, e);
        theta3 = ComputeTrueAnomalyFromMeanAnomaly(Me3, e);


        r1_ = a * ( 1 - e^2 ) / ( 1 + e* cos(theta1) );
        r2_ = a * ( 1 - e^2 ) / ( 1 + e* cos(theta2) );
        r3_ = a * ( 1 - e^2 ) / ( 1 + e* cos(theta3) );

        % 2) Compute Lagrance Coefficients

        f1 = ( 1 - mu * r1_/h^2 * ( 1 - cos(theta1 - theta2) ) + f1)/2;
        f3 = ( 1 - mu * r3_/h^2 * ( 1 - cos(theta3 - theta2) ) + f3)/2;
        g1 = ( r1_*r2_/h*sin(theta1-theta2) + g1)/2;
        g3 = ( r3_*r2_/h*sin(theta3-theta2) + g3)/2;

        % 3) Calculate c1, c3
        c1 = g3 / ( f1*g3 - f3*g1 );
        c3 = -g1 / ( f1*g3 - f3*g1 );

        % 4) Compute rho magnitudes
        rhoMag1 = ( -D(1,1) + D(2,1)/c1 - c3*D(3,1)/c1 ) / D0;
        rhoMag2 = ( -c1*D(1,2) + D(2,2) - c3*D(3,2) ) / D0;
        rhoMag3 = ( -c1*D(1,3)/c3 + D(2,3)/c3 - D(3,3) ) / D0;

        % 5) Compute r1, r2, r3
        r1 = RSite(:, 1) + rhoMag1 * rho(:, 1);
        r2 = RSite(:, 2) + rhoMag2 * rho(:, 2);
        r3 = RSite(:, 3) + rhoMag3 * rho(:, 3);

        % Compute absolute relative err
        rhos = [rhoMag1;rhoMag2;rhoMag3];
        
        % Compute relative error on each element
        relErrVec = (rhos./oldRhos - 1);
        
        % Compute magnitude of relative error vector
        absRelativeError = sqrt(relErrVec'*relErrVec);
        
        % Update for next iteration
        oldRhos = [rhoMag1;rhoMag2;rhoMag3];
        
        % Increment number of iterations
        nIter = nIter + 1;
    end
    
end

% 11) Apply Gibs Method
[a, e, i, RAAN, w, TA, h, Eps, T] = GibbsOrbitDetermination(r1, r2, r3, mu);





