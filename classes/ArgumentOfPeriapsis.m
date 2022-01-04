% This class contains methods for computing ArgumentOfPeriapsis
classdef ArgumentOfPeriapsis
    
    
    %% Generic Functions
    methods (Static)
        % Given specific angular momentum vector and eccentricity vector in Eci
        function ome = ComputeFrom_hVecEci_eVecEci(hVecEci, eVecEci)
            % Eci z is defined as
            zEci = [0;0;1];
            
            n = vec_crossp(zEci, hVecEci);
            
            % Compute ome
            ome = acos( dot(n,eVecEci) ./ (vec_mag(n) .* vec_mag(eVecEci)) );
            
            % Complete to 360deg(unwrap)
            ome((eVecEci'*zEci) < 0) = 2*pi -  ome((eVecEci'*zEci) < 0);
        end
        
    end
    %% Elliptic Orbit functions
    %% Parabolic Orbit functions
    %% Hyperbolic Orbit functions
end