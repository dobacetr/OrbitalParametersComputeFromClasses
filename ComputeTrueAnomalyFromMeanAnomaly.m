function trueAnomaly = ComputeTrueAnomalyFromMeanAnomaly(meanAnomaly, eccentricity)

%Me = E - e sinE
%Me - E + e * sinE = 0 = f
%df / dE = -1 + e * cosE

% Maximum number of iterations
iterMax = 100;

% Initial Guess
eccentricAnomaly = meanAnomaly;

% Compute value
f = meanAnomaly - eccentricAnomaly + eccentricity * sin(eccentricAnomaly);

% Threshold
threshold = 1E-5 * pi/180;

% number of iterations
iterNo = 0;

while abs(f) > threshold && iterNo < iterMax
    % Increment number of iterations performed
    iterNo = iterNo + 1;
    
    % Compute Derivative
    dfdE = -1 + eccentricity * cos(eccentricAnomaly);

    % Compute step to f = 0
    step = -f/dfdE;

    % Take a step
    eccentricAnomaly = eccentricAnomaly + step;

    % Compute new value
    f = meanAnomaly - eccentricAnomaly + eccentricity * sin(eccentricAnomaly);

end

trueAnomaly = 2*atan(tan(eccentricAnomaly/2) * sqrt( (1+eccentricity)/(1-eccentricity) ));

end
