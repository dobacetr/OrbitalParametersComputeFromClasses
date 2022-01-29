function psi = Compute_Residuals_Jn(R,Jn,rhoMag,c,d)
% Convert magnitude vector into column vector
rhoMag = reshape(rhoMag, [], 1);

% Initialize residuals
psi = Jn*rhoMag;

% Add the remaining terms
for i = 1 : size(Jn, 2)-2
    psi( 3 * ( i-1 )+1 : 3 * i, 1 ) = psi( 3 * ( i-1 )+1 : 3 * i, 1 ) + c(i) * R(:, i) - R(:, i+1) + d(i) * R(:, i+2);
end
