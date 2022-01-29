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
    
    %% Hyperbolic Functions
    methods (Static)
        % Given semi-major radius(a), eccentricity(e), hyperbolic mean
        % anomaly(F)
        function r = ComputeFrom_a_e_F(a, e, F)
            assert( all(e>1), 'This equation is valid for hyperbolic orbits only.' )
            
            r = a .* ( 1 - e .* cosh(F) );
            
        end
    end
end