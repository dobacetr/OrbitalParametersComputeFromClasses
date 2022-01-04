% This class contains methods for computing MeanAnomaly
classdef MeanAnomaly
    
    
    %% Generic Functions
    %% Elliptic Orbit functions
    methods (Static)
        % Given Eccentric Anomaly(E), eccentricity(e)
        function Me = ComputeFrom_E_e(E, e)
            % This function is only valid for elliptic orbits
            assert(all(e<1), 'This function is valid for only elliptic orbits.');
            
            Me = E - e * sin(E);
        end
        
        % Given Period(T), time elapsed since periapsis(t)
        function Me = ComputeFrom_T_t(T, t)
            Me = 2*pi*t./T;
        end
    end
    %% Parabolic Orbit functions
    methods (Static)
        % Given True Anomaly(ta)
        function Mp = ComputeFrom_ta(ta, e)
            % This function is only valid for parabolic orbit
            assert( all(e==1), 'This function is valid for only parabolic orbits.' );
            
            % Barker's Equation
            Mp = 0.5 * tan(ta/2) + tan(ta/2).^3/6;
        end
        
        % Given time from periapsis(t), gravitational parameter(mu),
        % specific angular momentum(h)
        function Mp = ComputeFrom_t_mu_h(t, mu, h)
            Mp = mu.^2./h.^3.*t;
        end
    end
    %% Hyperbolic Orbit functions
    methods (Static)
        % Given Eccentric Anomaly(F), eccentricity(e)
        function Mh = ComputeFrom_F_e(F, e)
            Mh = e * sinh(F) - F;
        end
        
        % Given time from periapsis(t), specific angular momentum(h),
        % gravitational parameter(mu), eccentricity(e)
        function Mh = ComputeFrom_t_h_mu_e(t, h, mu, e)
            % This function is only valid for hyperbolic orbits
            assert( all(e>1), 'This function is valid for only hyperbolic orbits.' );
            
            Mh = mu.^2 ./ h.^3 .* t .* ( e.^2 - 1 ) .^(3/2);
        end
    end
end