function rho = ComputeObservationDirection(RA, DEC)
% Computes the observation vector direction from the site of measurement
cDEC = cos(DEC);
sDEC = sin(DEC);
cRA  = cos(RA);
sRA  = sin(RA);
rho = [cRA.*cDEC; cDEC.*sRA; sDEC];
end