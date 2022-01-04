% This class contains methods for computing Inclination
classdef Inclination
    
    
    %% Generic Functions
    methods (Static)
        % Given specific angular momentum vector in Eci
        function i = ComputeFrom_hVecEci(hVecEci)
            % zEci is defined as
            zEci = [0;0;1];
            
            % Compute i
            i = acos( (hVecEci' * zEci)' ./ vec_mag(hVecEci) );
        end
        
    end
    %% Elliptic Orbit functions
    %% Parabolic Orbit functions
    %% Hyperbolic Orbit functions
end