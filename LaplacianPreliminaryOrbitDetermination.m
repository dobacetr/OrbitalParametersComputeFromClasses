function [a, e, i, RAAN, w, TA, h, Eps, T] = LaplacianPreliminaryOrbitDetermination(JD, time, RA, DEC, idx, dIdx)

RA = reshape(RA, 1, []);
DEC = reshape(DEC, 1, []);

% Earth [m^3/s^2]
muEarth = 0.3986004418E6 * (1E3)^3;
muSun = 132712E6*(1E3)^3;

% Reference julian date 1 Jan 2000 12:00:00
JDt0 = 2459545.00000;

% Julian years from reference date
T = (JD - JDt0)/365.25;

% Mean Longitude of Earth at J2000
Me0 = 100.46435 * pi/180;

% Mean Longitude at time JD
Me = Me0 + T * 2*pi;

% Earts orbit eccentricity
e = 0.01671022;

% True Longitude of Earth at time JD
Lon = ComputeTrueAnomalyFromMeanAnomaly(Me, e);

% Earths Orbits Semi-Major Axis
a = 149597870700 * 1.00000011;

% Earths distance from sun
r = a * ( 1 - e^2 ) / ( 1 + e* cos(Lon) );

% Earths position in Barycentric Inertial Frame
C_perifocal2Bary = dcm_eci2perifocal(-11.26064*pi/180, 0.00005*pi/180, 102.94719*pi/180)';

% Earths position wrt sun in perifocal frame
rP = r*[cos(Lon); sin(Lon);0];

% Earths position in barycentric frame
q = C_perifocal2Bary * rP;

% Earths position wrt sun
qMag = sqrt(q'*q);
qDir = q/qMag;

% Define proper motion
dxdt = @(x, idx, dIdx) diff(x(:, idx+[-1,1]*dIdx), 1, 2) / diff( time(idx+[-1,1]*dIdx) );
% RAdt    = dxdt(RA, idx, dIdx);
% DECdt   = dxdt(DEC, idx, dIdx);
% dsdt = sqrt( RAdt * cos(DEC(idx))^2 + DECdt^2 );
RAdtAt    = @(idx, dIdx) dxdt(RA, idx, dIdx);
DECdtAt   = @(idx, dIdx) dxdt(DEC, idx, dIdx);
dsdtAt = @(idx, dIdx) sqrt( RAdtAt(idx,dIdx) * cos(DEC(idx))^2 + DECdtAt(idx, dIdx)^2 );
dsdt = dsdtAt(idx, dIdx);

dxds = @(x, idx, dIdx) dxdt(x, idx, dIdx) / dsdt;

% Moving orthonormal basis
rhoDir = ComputeObservationDirection(RA(idx), DEC(idx));

% vDir = derho/ds = derho/dt*dt/ds;
vDirAtIdx = @(idx, dIdx) dxds(ComputeObservationDirection(RA, DEC), idx, dIdx);
vDir = vDirAtIdx(idx, dIdx);
nDir = cross(rhoDir, vDir);

% Compute Kappa, geodesic curvature
% dvDir/ds = -rhoDir + K * nDir
% dvDir/ds = dvDir/dt * dt/ds
dvDirds = diff([vDirAtIdx(idx-dIdx, dIdx), vDirAtIdx(idx+dIdx, dIdx)], 1, 2)/diff(time(idx+[-1,1]*dIdx)) / dsdt;
K = (dvDirds + rhoDir) ./ nDir;
error('K can not be computed this way.');

% Compute C
C = dsdt^2*K*qMag^3/( muSun * dot(qDir, nDir) );

% Compute cosEps
cEps = dot( qDir, rhoDir );

% 8th order eq
fun = @(r) C^2 * r.^8 - qMag^2*(C^2 + 2*C*cEps+1)*r.^6+2*qMag^5*(C*cEps+1)*r.^3-q^8; % = 0

initGuess = r; % At Earth

% Try to find a real root
errStatus = 1;
r2 = 0;
% We dont check for real r2 because NewtonRaphson can not compute complex
% number
while ( ...
        errStatus ~= 0 ...      % Check for errors
      ) ...
        && ...
      ( ...
        initGuess < r+100E6 ...   % Check for initial guess attemps 
      )
    
%     [r2, residual, errStatus] = NewtonRaphsonRoot(fun, initGuess, 1E-8, 1.0, 100, @(x) x < (r-200E6));
    [r2, residual, errStatus] = SecondOrderNewtonRaphsonRoot(fun, initGuess, 1E-8, 1, 100, @(x) x < (r-200E6));
    initGuess = initGuess + 1E6;
end
% roots([1, 0, a, 0, 0, b, 0 , 0, c])
assert( initGuess < r+100E6, 'Exceeded limit for initial guess to NewtonRapson.' );

% compute rhoMag
rhoMag = (1-(qMag/r2)^3)*q/C;

d2sdt2 = diff([dsdtAt(idx-dIdx, dIdx), dsdtAt(idx+dIdx, dIdx)], 1 ,2) / time(idx+[-1,1]*dIdx);
rhoMagDot =  (muSun * dot( q, vDir ) * ( 1/qMag^3 - 1/r2^3 ) - rhoMag*d2sdt2)/(2*dsdt);

v2 = rhoMag * vDir + rhoMagDot * rhoDir;

[a, e, i, RAAN, w, TA, h, Eps, T] = ComputeKepplerianElements(r2, v2, muEarth);
end