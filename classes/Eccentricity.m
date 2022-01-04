% This is a class that containts methods for computing eccentricity
classdef Eccentricity
    
    %% Generic Orbit functions
    methods (Static)
                
        % Given position(Rvec) and velocity(Vvec) as vectors and gravitational
        % parameter(Mu)
        function eVec = ComputeFrom_rVec_vVec_Mu(rVec, vVec, mu)
            
            % Compute specific angular momentum
            hVec = SpecificAngularMomentum.ComputeFrom_rVec_vVec(rVec, vVec);
            
            % Compute C
            C = cross(vVec, hVec) - mu .* rVec ./ vec_mag(rVec);
            
            % Compute eccentricty vector
            eVec = C ./ mu;
            
        end
        
        % Given range(r), specific angular momentum(h), true anomaly(ta) and
        % gravitational parameter(mu)
        function e = ComputeFrom_r_h_ta_mu(r, h, ta, mu)
            % r = h^2/mu * ( 1 / ( 1 + e*cos(ta) ) )
            % Solve for eccentricity
            e = ( h.^2./(mu.*r) - 1 ) ./ cos(ta);
            
            % Verify e is positive
            if isa(e, 'sym')
            else
                assert( all(e>=0), 'Eccentricity must be positive.' );
            end
        end
                
        % Given periapsis range(rp), specific angular momentum(h),
        % gravitational parameter(mu)
        function e = ComputeFrom_rp_h_mu(rp, h, mu)
            % at rp, we have ta = 0;
            e = Eccentricity.ComputeFrom_r_h_ta_mu(rp, h, 0, mu);
        end
        
    end
    %% Elliptic Orbit functions
    methods (Static)

        % Given periapsis and apoapsis ranges
        function e = ComputeFrom_rp_ra(rp, ra)
            % This is valid from elliptical orbits only. Parabolic and
            % hyperbolic orbits should have infinity as their ra.
            
            % Check that ra > rp
            if isa(ra,'sym') || isa(rp,'sym')
            else
                assert( all( ra >= rp ), ' Apoapsis must be greater than periapsis.' );
            end
            
            % Compute eccentricity
            e = ( ra - rp ) ./ ( ra + rp );
            
            % Verify that the orbit is elliptical
            if isa(e, 'sym')
            else
                assert( all( e <1 ), 'This function is valid for elliptical orbits only.' );
            end
        end

        % Given apoapsis range(rp), specific angular momentum(h),
        % gravitational parameter(mu)
        function e = ComputeFrom_ra_h_mu(ra, h, mu)
            % at rp, we have ta = 180deg;
            e = ComputeFrom_r_h_ta_mu(ra, h, pi, mu);
        end
        
        % Given periapsis range(rp), semi-major axis(a)
        function e = ComputeFrom_rp_a(rp, a)
            % This is only valid from elliptical orbits
            e = 1 - rp./a;
            
            % Check if the orbits are elliptical
            assert( all( e < 1 ), 'This function is valid for elliptical orbits only.' );
        end
        
        % Given apoapsis range(rp), semi-major axis(a)
        function e = ComputeFrom_ra_a(ra, a)
            % This is only valid for elliptical orbits.
            % Compute periapsis
            rp = 2 * a - ra;
            
            % Compute from periapsis and semi-major axis
            e = ComputeFrom_rp_a(rp, a);
        end
        
        % Given semi-major axis(a), specific angular momentum(h),
        % gravitational parameter(mu)
        function e = ComputeFrom_a_h_mu(a, h, mu)
            % This function is only valid for elliptical orbits
            e2 = 1 - h.^2 ./ ( mu .* a );
            
            % Check for real valued eccentricity and elliptical orbit
            assert( all(e2 > 0) && all( e2 < 1 ), 'This function is valid for elliptical orbits only.' );
            
            % Compute eccentricity
            e = sqrt(e2);
        end
        
    end
    %% Parabolic Orbit functions
    %% Hyperbolic Orbit functions
end