function [c, d] = Compute_c_d_Jn(fkm1, fkp1, gkm1, gkp1)
% Computes coefficients for the Jacobian of Jn Method
c =  gkp1 ./ ( fkm1.*gkp1 - fkp1.*gkm1 );
d = -gkm1 ./ ( fkm1.*gkp1 - fkp1.*gkm1 );
end