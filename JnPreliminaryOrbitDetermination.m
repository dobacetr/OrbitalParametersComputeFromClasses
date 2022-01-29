function  [a, e, i, RAAN, w, TA, h, Eps, T] = JnPreliminaryOrbitDetermination(time, rho, RSiteECI, mu)

% For simplicity, assume constant range for all measurements
rhoInit = 5000E3;

% Initialize absolute relative error
absErrDiff = Inf;
absErr    = Inf;
absErrOld = Inf;

% Set stopping tolerance, residual and maximum iterations
tolAbsErrDiff = 1E-5;
tolResidual = 5E-4;
nIterMax = 1E5;

% Initialize iteration counter
nIter = 0;

% Check time difference between data. f, g computations were written with
% constant time intervals.
assert(all(abs(diff(time,2))<1E-10), 'The f, g methods was written with equitime interval condition. Check time values');
dt = time(2) - time(1);
    
% Initial range estimates
rhoMag = ones( size(rho, 2), 1 ) * rhoInit;

% Weight matrix
W = eye(3 * (size(rho, 2)- 2));

% step scale
scale = 0.5;

while nIter < nIterMax ...
        && ( ...
                absErrDiff > tolAbsErrDiff ...
                || absErr > tolResidual ...
                )
        

    % Compute position vectors with the current estimated ranges
    r = RSiteECI + rho .* rhoMag';
    
    % Compute Velocity trough lower order lagrange coefficients
    [fkm1, fkp1, gkm1, gkp1] = Compute_f_g_from_r_dt(r(:, 2:end-1), dt, mu);
    rkm1 = r(:, 1:end-2);
    rkp1 = r(:,3:end);
    v = -fkp1./(fkm1.*gkp1-fkp1.*gkm1).*rkm1 + fkm1./(fkm1.*gkp1-fkp1.*gkm1).*rkp1;
    v = [v(:, 1), v, v(:, end)]; % Pad the beginning and ending

    % Compute Jacobian Matrix trough higher order lagrange coefficients
    [fkm1, fkp1, gkm1, gkp1] = Compute_f_g_from_r_v_dt(r, v, dt, mu);
    [c, d] = Compute_c_d_Jn(fkm1, fkp1, gkm1, gkp1);
    Jn = ConstructJacobian_Jn(c,d,rho);
    
    % Compute Residuals
    psi = Compute_Residuals_Jn(RSiteECI,Jn,rhoMag,c,d);

    % Compute required change in ranges
    dRhoMag = - (Jn'*W*Jn)\Jn'*W*psi;
    
    % Compute new ranges
    rhoMag = rhoMag + dRhoMag*scale;
    
    % Compute Residual Magnitude
    absErr = sqrt(psi'*psi);
    
    % Compute Residual Magnitude Change
    absErrDiff = abs(absErr-absErrOld);
    
    % Increment iteration
    nIter = nIter + 1;
    
    % Save absErr for next iteration
    absErrOld = absErr;
    
%   Record for debugging
%     hist_absErr(nIter) = absErr;
%     hist_absErrDif(nIter) = absErrDiff;
%     hist_rhoMag(:, nIter) = rhoMag;
%     hist_dRhoMag(:, nIter) = dRhoMag;
end

if absErr > tolResidual
    error('Solution did not converge.')
end

% Compute position vectors
r = RSiteECI + rho .* rhoMag';

% TODO: employ a mean orbit estimation that uses all positions
% Compute orbital parameters from the first 3 positions
[a, e, i, RAAN, w, TA, h, Eps, T] = GibbsOrbitDetermination(r(:, 1), r(:, 2), r(:, 3), mu);
end