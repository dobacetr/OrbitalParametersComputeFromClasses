% This class contains methods for computing EccentricAnomaly
classdef EccentricAnomaly
    
    
    %% Generic Functions
    methods (Static)
    end
    %% Elliptic Orbit Functions
    methods (Static)
        % Given MeanAnomaly(Me) and eccentricity(e)
        function E = ComputeFrom_Me_e(Me, e)
            assert(all(e<1), 'This function is valid for only elliptic orbits.');
            
            assert( ...
                length(Me) == 1 || ...
                length(e) == 1 || ...
                length(Me) == length(e), ...
                'Number of elements of Me and e must match.');
            
            sz = max(length(Me), length(e));
            
            % Initialize results array
            E = zeros(1, sz);
            
            % Indices for Me and e if they are vectors
            idxMe = 1;
            idxe  = 1;
            
            % Solve for each point independently
            for iRes = 1 : sz
                
                % Determine index of Me to be used
                if length(Me) > 1
                    idxMe = iRes;
                end
                
                % Determine index of e to be used
                if length(e) > 1
                    idxe = iRes;
                end
                
                % Function to be solved
                %Me = E - e sinE
                %Me - E + e * sinE = 0 = f
                fun = @(E) Me(idxMe) - E + e(idxe) * sin(E);
                
                % Initial Guess
                Einit = Me(idxMe);
                
                [E(iRes), residual, errStatus] = SecondOrderNewtonRaphsonRoot(fun, Einit);
                
                if errStatus > 0
                    error('Error occurred while searching for Eccentric Anomaly.');
                end
            end
        end
        
        
        % Given true anomaly(ta), eccentricity(e)
        function E = ComputeEFrom_ta_e(ta, e)
            E = 2 * atan( sqrt( (1-e)./(1+e) ) .* tan(ta/2) );
        end
    end
    
    %% Parabolic Orbit Functions
    %% Hyperbolic Orbit Functions
    methods (Static)
        % Given MeanAnomaly(Mh) and eccentricity(e)
        function F = ComputeFrom_Mh_e(Mh, e)
            assert(all(e>1), 'This function is valid for only hyperbolic orbits.');
            
            assert( ...
                length(Mh) == 1 || ...
                length(e) == 1 || ...
                length(Mh) == length(e), ...
                'Number of elements of Mh and e must match.');
            
            sz = max(length(Mh), length(e));
            
            % Initialize results array
            F = zeros(1, sz);
            
            % Indices for Mh and e if they are vectors
            idxMh = 1;
            idxe  = 1;
            
            % Solve for each point independently
            for iRes = 1 : sz
                
                % Determine index of Me to be used
                if length(Mh) > 1
                    idxMh = iRes;
                end
                
                % Determine index of e to be used
                if length(e) > 1
                    idxe = iRes;
                end
                
                % Function to be solved
                %Mh = e sinF - F
                %Mh - e * sinF + F = 0 = f
                fun = @(F) Mh(idxMh) - e(idxe) * sinh(F) + F;
                
                % Initial Guess
                Einit = Mh(idxMh);
                
                [F(iRes), residual, errStatus] = SecondOrderNewtonRaphsonRoot(fun, Einit);
                
                if errStatus > 0
                    error('Error occurred while searching for Eccentric Anomaly.');
                end
            end
        end
        
        % Given true anomaly(ta), eccentricity(e)
        function F = ComputeFFrom_ta_e(ta, e)
            F = 2 * atanh( sqrt( (e-1)./(1+e) ) .* tan(ta/2) );
        end
    end
end