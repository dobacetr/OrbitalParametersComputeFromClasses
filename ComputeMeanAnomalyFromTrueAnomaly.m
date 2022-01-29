function meanAnomaly = ComputeMeanAnomalyFromTrueAnomaly(trueAnomaly, eccentricity)
% Compute Eccentric Anomaly
E = 2 * atan( sqrt( (1-eccentricity)/(1+eccentricity) )*tan(trueAnomaly/2) );

% Compute Mean Anomaly
meanAnomaly = E - eccentricity * sin(E);