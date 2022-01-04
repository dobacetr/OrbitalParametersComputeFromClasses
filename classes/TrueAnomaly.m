% This class contains methods for computing TrueAnomaly
classdef TrueAnomaly
    
    
    %% Generic Functions
    methods (Static)
        % Given specific angular momentum vector and eccentricity vector in Eci
        function ta = ComputeFrom_rVecEci_eVecEci(rVecEci, eVecEci)
            % Eci z is defined as
            zEci = [0;0;1];
            
            % Compute ta
            ta = acos( dot(rVecEci,eVecEci) ./ (vec_mag(rVecEci) .* vec_mag(eVecEci)) );
            
            % Complete to 360deg(unwrap)
            ta((rVecEci'*zEci) < 0) = 2*pi -  ta((rVecEci'*zEci) < 0);
        end
        
        % Given range(r), specific angular momentum(h), eccentricity(e),
        % gravitational parameter(mu) and 
        % whether 0 <= ta <= pi or pi < ta < 2*pi
        function ta = ComputeFrom_r_h_e_mu(r, h, e, mu, isPastApoapsis)
            ta = acos( ( h.^2 ./ ( mu .* r) -1 ) ./ e );
            
            ta(isPastApoapsis) = 2*pi - ta(isPastApoapsis);
        end
    end
    %% Elliptic Orbit functions
    methods (Static)
        % Given range(r), semi-major axis(a), eccentricity(e), and whether
        % 0 <= ta <= pi or pi < ta < 2*pi
        function ta = ComputeFrom_r_a_e(r,a,e, isPastApoapsis)
            ta = acos( ( a .* ( 1 - e.^2 ) ./ r - 1 ) ./ e );
            
            ta(isPastApoapsis) = 2*pi - ta(isPastApoapsis);
        end
        
        % Given MeanAnomaly(Me), eccentricty(e)
        function ta = ComputeFrom_Me_e(Me, e)
            % Compute Eccentric Anomaly
            E = EccentricAnomaly.ComputeFrom_Me_e(Me, e);
            
            % Compute True Anomaly
            ta = ComputeFrom_E_e(E, e);
        end
        
        % Given EccentricAnomaly(E), eccentricity
        function ta = ComputeFrom_E_e(E, e)
            ta = 2*atan(tan(E/2) .* sqrt( (1+e)./(1-e) ));
        end
    end
    %% Parabolic Orbit functions
    methods (Static)
        % Given MeanAnomaly(Mp)
        function ta = ComputeFrom_Mp(Mp, e)
            % This function is only valid for parabolic orbits
            assert( all(e==1), 'This function is valid for only parabolic orbits.' );
            
            temp = ( 3*Mp + sqrt( (3*Mp).^2 + 1 ) );
            
            % Use the sign Mp to determine the sign of ta
            ta = sign(Mp) .* 2 .* atan( temp^(1/3) + temp^(-1/3)  );
        end
    end
    %% Hyperbolic Orbit functions
    methods (Static)
        % Given MeanAnomaly(Mh), eccentricity(e)
        function ta = ComputeFrom_Mh_e(Mh, e)
            % Compute Eccentric Anomaly
            F = EccentricAnomaly.ComputeFrom_Mh_e(Mh, e);
            
            % Use the sign Mh to determine the sign of ta
            ta = sign(Mh) .* 2 .* atan( tanh(F/2) .* sqrt( (e+1)./(e-1) ) );
        end
    end
end