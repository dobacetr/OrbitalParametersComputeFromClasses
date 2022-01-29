function [a, e, i, RAAN, w, TA, h, Eps, T] = ComputeKepplerianElements(ri, vi, mu)
% This function computes kepplerian elements
% It should be noted at first 6 of these are enough to define the orbit
% completely.

% ri, vi are in inertial frame

% X of inertial frame
xi = [1;0;0];

% Y of inertial frame
yi = [0;1;0];

% Z of inertial frame
zi = [0;0;1];

% Specific angular momentum
hVec = cross(ri, vi);
h = sqrt(hVec' * hVec);

% eccentricity
rMag = sqrt(ri' * ri);
eVec = cross(vi, hVec)/mu - ri / rMag;
e = sqrt(eVec' * eVec);

% semi-major axis
a = h^2/(mu * ( 1 - e^2 ) );

% inclination
i = acos(  dot(hVec, zi) / h );

% specific energy
Eps = dot(vi, vi) /2 - mu / rMag;

% Right ascension of the ascending one
n = cross( zi, hVec );
nMag = sqrt(n'* n);
RAAN = acos( dot( n, xi ) / nMag );

if dot(n, yi) >= 0
    % Do Nothing, 0<=RAAN<=pi
else
    RAAN = 2*pi - RAAN;
end

% argument of perigee
w = acos( dot(n, eVec)/ ( nMag * e ) );
if dot(eVec, zi) >= 0
    % Do Nothing, 0<= w <= pi
else
    w = 2*pi - w;
end

% true anomaly
TA = acos( dot(ri, eVec) / ( rMag * e ) );
if dot(ri, zi) >= 0
    % Do Nothing, 0<= TA <= pi
else
    TA = 2*pi - TA;
end

% period
T = 2*pi/sqrt(mu)*a^(1.5);

end