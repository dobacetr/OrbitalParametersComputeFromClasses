% This class contains methods for computing SpecificEnergy
classdef SpecificEnergy
    
    
    %% Generic Functions
    methods (Static)
        % Given position vector(rVec), velocity vector(vVec), mu
        function eps = ComputeFrom_rVec_vVec_mu(rVec, vVec, mu)
            eps = dot(vVec, vVec) / 2 - mu ./ vec_mag(rVec);
        end
        
        % Given position vector(rVec), velocity vector(vVec), mu
        function eps = ComputeFrom_r_v_mu(r, v, mu)
            eps = v.^2 / 2 - mu ./ r;
        end
        
        % Given semi-major axis(a), mu
        function eps = ComputeFrom_a_mu(a, mu, e)
            % If the orbit type is known, a random e within its range can
            % be given as it is not used in the computation
            if a < 0
                % Semi-major axis is negative, then e > 1
                assert( e > 1, 'Semi-major axis is negative while orbit is not hyperbolic.' );

            elseif e > 1
                % if a is actually abs(a), but the orbit is hyperbolic
                a = -a;
                
            else
                % Orbit is not hyperbolic, there should not be any
                % confusion about semi-major axis
                %
                % Do nothing
            end
                
            eps = -mu ./ (2*a);
        end
        
        % Given specific angular momentum(h), eccentricity(e),
        % gravitational parameter(mu)
        function eps = ComputeFrom_h_e_mu(h, e, mu)
            eps = -0.5 * mu.^2 ./ h.^2 * ( 1 - e.^2 );
        end
    end
    %% Elliptic Orbit functions
    %% Parabolic Orbit functions
    %% Hyperbolic Orbit functions
end