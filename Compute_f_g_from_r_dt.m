function [fkm1, fkp1, gkm1, gkp1] = Compute_f_g_from_r_dt(r, dt, mu)
% Computes 4th order fk-1 fk+1 gk-1 gk+1 from equidistance points, using rk
% vk
r_ = sqrt(sum(r.*r, 1));
fkm1 = 1 - mu./(2*r_.^3)*dt^2 ;
fkp1 = 1 - mu./(2*r_.^3)*dt^2 ;
gkm1 = -dt + mu./(6*r_.^3)*dt^3;
gkp1 = dt - mu./(6*r_.^3)*dt^3;
end