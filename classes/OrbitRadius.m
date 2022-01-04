% This class contains methods for computing OrbitRadius
classdef OrbitRadius
    
    
    %% Generic Functions
    methods (Static)
        % Given specific angular momentum(h), eccentricity(e), true
        % anomaly(ta), gravitational parameter(mu)
        function r = ComputeFrom_h_e_ta_mu(h, e, ta, mu)
            r = h.^2 ./ ( mu .* ( 1 + e .* cos(ta) ) );
        end
    end
end