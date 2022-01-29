function [fkm1, fkp1, gkm1, gkp1] = Compute_f_g_from_r_v_dt(r, v, dt, mu)
% Computes 4th order fk-1 fk+1 gk-1 gk+1 from equidistance points, using rk
% vk
r_ = sqrt(sum(r.*r, 1));
v_ = sqrt(sum(v.*v, 1));
fkm1 = 1 - mu./(2*r_.^3)*dt^2 - mu * dot(r, v) ./ ( 2 * r_ .^5 ) * dt^3 + mu/24*( -2*mu./r_.^6+3*v_.^2./r_.^5 - 15 * dot(r,v).^2./r_.^7 )*dt^4;
fkp1 = 1 - mu./(2*r_.^3)*dt^2 + mu * dot(r, v) ./ ( 2 * r_ .^5 ) * dt^3 + mu/24*( -2*mu./r_.^6+3*v_.^2./r_.^5 - 15 * dot(r,v).^2./r_.^7 )*dt^4;
gkm1 = -dt + mu./(6*r_.^3)*dt^3+mu*dot(r, v)./(4*r_.^5)*dt^4;
gkp1 = dt - mu./(6*r_.^3)*dt^3+mu*dot(r, v)./(4*r_.^5)*dt^4;
end