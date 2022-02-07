function [a, e, i, RAAN, w, TA, h, Eps, T, v2] = GibbsOrbitDetermination(r1, r2, r3, mu)
% It should be noted at first 6 of these are enough to define the orbit
% completely.

% 1) r1, r2, r3 are given

% 2) Compute pairwise planes
C12 = cross(r1, r2);
C23 = cross(r2, r3);
C31 = cross(r3, r1);

% Tolerance for coplanarity
tol = 1E-6;
coPlanarity = abs(dot(r1, C23)) / sqrt( (r1'*r2) * (C23'*C23) );

assert( coPlanarity < tol, 'Position vectors are coplanar.' );

% 3)
r1Mag = sqrt(r1'*r1);
r2Mag = sqrt(r2'*r2);
r3Mag = sqrt(r3'*r3);

N = r1Mag * C23 + r2Mag*C31 + r3Mag * C12;
D = C12 + C23 + C31;
S = r1 * (r2Mag-r3Mag) + r2 * ( r3Mag-r1Mag ) + r3 * (r1Mag-r2Mag);

% 4) Compute V2
NMag = sqrt( N'*N );
DMag = sqrt( D'*D );

v2 = sqrt( mu / ( NMag*DMag ) ) * ( cross(D, r2)/r2Mag + S );

% Compute and return Kepplerian orbit elements
[a, e, i, RAAN, w, TA, h, Eps, T] = ComputeKepplerianElements(r2, v2, mu);
end