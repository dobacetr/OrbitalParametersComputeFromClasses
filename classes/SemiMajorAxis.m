% This class contains methods for computing Semi Major Axis
classdef SemiMajorAxis
    
    
    %% Generic Functions
    methods (Static)
        % Given specific energy(eps), gravitational parameter(mu)
        function a = ComputeFrom_eps_mu(eps, mu)
            a = - mu ./ ( 2*eps );
        end
        
        % Given position vector(rVec), velocity vector(vVec), gravitational
        % parameter(mu)
        function a = ComputeFrom_rVec_vVec_mu(rVec, vVec, mu)
            % Compute specific energy
            eps = SpecificEnergy.ComputeFrom_rVec_vVec_mu(rVec, vVec, mu);
            
            % Compute a
            a = SemiMajorAxis.ComputeFrom_eps_mu(eps, mu);
        end
        
        % Given specific angular momentum(h), eccentricity(e),
        % gravitational parameter(mu)
        function a = ComputeFrom_h_e_mu(h, e, mu)
            a = h.^2./(mu .* ( 1 - e.^2 ) );
        end
    end
    %% Elliptic Orbit functions
    methods (Static)
        
        % Given periapsis range(rp), apoapsis range(ra)
        function a = ComputeFrom_rp_ra(rp, ra, e)
            % Note that you can give a random e<1, if you are sure it is an
            % elliptical orbit, without actually knowing the value.
            
            % This is only valid for elliptical orbits
            if isa(e, 'sym')
            else
                assert( all(e < 1), 'This function is valid for only elliptical orbits' );
            end
            
            % Compute a
            a = ( rp + ra ) / 2;
        end

        % Given Period(T), gravitational parameter(mu)
        function a = ComputeFrom_T_mu(T, mu, e)
            % Note that you can give a random e<1, if you are sure it is an
            % elliptical orbit, without actually knowing the value.
            
            % This is only valid for elliptical orbits
            assert( all(e < 1), 'This function is valid for only elliptical orbits' );
            
            a = ( T / ( 2*pi ) .* sqrt(mu) ).^(2/3);
        end
    end
    %% Parabolic Orbit functions
    %% Hyperbolic Orbit functions
    methods (Static)
         
        % Given eccentricity(e), gravitational parameter(mu) for hyperbolic
        % orbits
        function a = ComputeFrom_mu_e(mu, e)
            % Check that the orbits are hyperbolic
            assert( all(e>1), 'This function is valid only for hyperbolic orbits.' );
            
            % Compute a
            a = mu ./ ( 2*e );
            
        end
    end
end