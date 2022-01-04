% This is class that contains methods for computing specific angular
% momentum
classdef SpecificAngularMomentum
    
    %% Generic Orbit functions
    methods (Static)
        
        % Given position and velocity as a vector
        function hVec = ComputeFrom_rVec_vVec(rVec, vVec)
            hVec = cross(rVec, vVec);
        end
        
        % Given position magnitude and tangential velocity
        function h = ComputeFrom_r_vTangent(r, vTangent)
            h = r.*vTangent;
        end
        
        
    end
    %% Elliptic Orbit functions
    methods (Static)
        
        % Given Semi-Major axis(a), eccentricity(e), gravitational
        % parameter(mu)
        function h = ComputeFrom_a_e_mu(a, e, mu)
            if isa(e, 'sym')
            else
                assert(all(e<1), 'This function is valid for only elliptic orbits.');
            end
            
            h = sqrt( a .* mu .* ( 1 - e.^2 ) );
        end
    end
    %% Parabolic Orbit functions
    %% Hyperbolic Orbit functions
    
end