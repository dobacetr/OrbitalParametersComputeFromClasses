function  [a, e, i, RAAN, w, TA, h, Eps, T] = JnPreliminaryOrbitDetermination(time, rho, RSiteECI, mu)

% For simplicity, assume constant range for all measurements
rhoInit = 5000E3;

% Initialize absolute relative error
absRelErr = Inf;
absErr    = Inf;

% Set stopping tolerance, residual and maximum iterations
tolAbsRelErr = 1;
tolResidual = 1E-3;
nIterMax = 1E5;
nAdaptStepSizeAt = 5;

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

% initial
oldDRhoMag = 1;
oldAbsRelErr = Inf;

% step scale
scale = 0.5;

% counter to decrease the step size
nCounterAdaptiveStep = 0;

while nIter < nIterMax ...
        && ( ...
                absRelErr > tolAbsRelErr ...
                || absErr > tolResidual ...
                )
        

    % Compute position vectors with the current estimated ranges
    r = RSiteECI + rho .* rhoMag';
    
    % Simple compute for velocity
    % Compute Lagrange coefficients to compute velocity
    [fkm1, fkp1, gkm1, gkp1] = Compute_f_g_from_r_dt(r(:, 2:end-1), dt, mu);
    rkm1 = r(:, 1:end-2);
    rkp1 = r(:,3:end);
    v = -fkp1./(fkm1.*gkp1-fkp1.*gkm1).*rkm1 + fkm1./(fkm1.*gkp1-fkp1.*gkm1).*rkp1;
    v = [v(:, 1), v, v(:, end)]; % Pad the beginning and ending
%     v = diff(r, 1, 2) ./ diff(time);
%     v = [v, v(:, end)];                 % Pad the end with the same velocity


    [fkm1, fkp1, gkm1, gkp1] = Compute_f_g_from_r_v_dt(r, v, dt, mu);
    [c, d] = Compute_c_d_Jn(fkm1, fkp1, gkm1, gkp1);
    Jn = ConstructJacobian_Jn(c,d,rho);
    psi = Compute_Residuals_Jn(RSiteECI,Jn,rhoMag,c,d);


    % Compute change in ranges
    dRhoMag = - (Jn'*W*Jn)\Jn'*W*psi;
    
    % Compute new ranges
    rhoMag = rhoMag + dRhoMag*scale;

    % Compute change in range vector
    relErrVec = (dRhoMag./oldDRhoMag-1);
    absRelErr = sqrt(relErrVec'*relErrVec);
    
    % Count consecutive non decreasing absRelErr
    if absRelErr > oldAbsRelErr &&...
            absErr < tolResidual % absErr converged
        nCounterAdaptiveStep = nCounterAdaptiveStep + 1;
    else
        nCounterAdaptiveStep = 0;
        % only save the lowest absRelErr
        oldAbsRelErr = absRelErr;
    end
    
    % When the counter reaches, reset
    if nCounterAdaptiveStep >= nAdaptStepSizeAt
        scale = scale * 0.5;
        nCounterAdaptiveStep = 0;
    end
    
    % Save the values for next iteration
    oldDRhoMag = dRhoMag;
    
    % Compute Residual
    absErr = sqrt(psi'*psi);
    
    % Increment iteration
    nIter = nIter + 1;
    
    hist_absErr(nIter) = absErr;
    hist_absRelErr(nIter) = absRelErr;
    hist_rhoMag(:, nIter) = rhoMag;
    hist_dRhoMag(:, nIter) = dRhoMag;
end

if absErr > tolResidual
    error('Solution did not converge.')
end

r = RSiteECI + rho .* rhoMag';
[a, e, i, RAAN, w, TA, h, Eps, T] = GibbsOrbitDetermination(r(:, 1), r(:, 2), r(:, 3), mu);
end