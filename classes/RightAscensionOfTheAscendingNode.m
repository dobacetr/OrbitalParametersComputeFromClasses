% This class contains methods for computing RightAscensionOfTheAscendingNode
classdef RightAscensionOfTheAscendingNode
    
    
    %% Generic Functions
    methods (Static)
        % Given specific angular momentum vector in Eci
        function Ome = ComputeFrom_hVecEci(hVecEci)
            % Eci x-y-z is defined as
            zEci = [0;0;1];
            yEci = [0;1;0];
            xEci = [1;0;0];
            
            n = vec_crossp(zEci, hVecEci);
            
            % Compute Ome
            Ome = acos( (n'*xEci)' ./ vec_mag(n) );
            
            % Complete to 360deg(unwrap)
            Ome((n'*yEci) < 0) = 2*pi -  Ome((n'*yEci) < 0);
        end
        
    end
    %% Elliptic Orbit functions
    %% Parabolic Orbit functions
    %% Hyperbolic Orbit functions
end