classdef OrbitVelocity
    %% Generic Orbit
    methods (Static)
        % Given radial velocity(vR), tangential velocity(vT), Right
        % Ascension of the Ascending Node(Ome), Inclination(i), Argument of
        % Perigee(w), True Anomaly(ta)
        % Compute eci velocity vector
        function vEci = ComputeFrom_vR_vT_Ome_i_w_ta(vR, vT, Ome, i, w, ta)
            
            % Velocity at RadialTangentialZ Coordinate System
            vRTZ = [vR;vT;zeros(1, length(vR))];
            
            % Compute dcm from perifocal to eci
            C_perifocal2eci = permute(dcm_eci2perifocal(Ome, i, w), [2,1,3]);
            
            % Compute dcm RTZ to perifocal
            sz = length(ta);
            cta = cos(ta);
            sta = sin(ta);
            C_RTZ2perifocal = reshape([...
                cta;
                sta;
                zeros(1, sz);
                -sta;
                cta;
                zeros(3, sz);
                ones(1, sz);
                ], 3,3, sz);
            
            % Compute Velocity at eci
            vEci = MVmult( Mmult(C_perifocal2eci , C_RTZ2perifocal), vRTZ);
        end
        
        function vEci = ComputeFrom_h_e_Ome_inc_w_ta_mu(h, e, Ome, inc, w, ta, mu)
            
            sze = length(e);
            szta = length(ta);
            
            assert( sze == szta || ( sze == 1 || szta == 1),'Incompatible eccentricity and true anomanly vector sizes.' );
            
            sz = max(sze, szta);
            
            vp = (mu./h) .* ([-sin(ta);(e + cos(ta));zeros(1, sz)]);
            
            Cip = permute(dcm_eci2perifocal(Ome, inc, w), [2,1,3]);
            
            vEci = MVmult(Cip, vp);
        end
        
        % Given specific energy(eps), radius(r), gravitational parameter(mu)
        function v = ComputeFrom_eps_r_mu(eps, r, mu)
            v = sqrt(2*( eps + mu ./ r ));
        end
    end
    %% Elliptic Orbit
    methods (Static)
        % given specific angular momentum(h), eccentricity(e), true
        % anomaly(ta) , gravitational parameter(mu)
        function [vTangential, vRadial] = ComputeFrom_h_e_ta_mu(h, e, ta, mu)
            % This might only be valid for elliptic orbits, so i will force
            % it
            if isa(e, 'sym')
            else
                assert( all(e<1), 'This function is valid for only elliptic orbits.' );
            end
            
            % Make sure 0<=ta<360deg
            ta = wrapTo2Pi(ta);
            
            % Compute radius
            r = OrbitRadius.ComputeFrom_h_e_ta_mu(h, e, ta, mu);
            
            % Compute Specific Energy
            eps = SpecificEnergy.ComputeFrom_h_e_mu(h, e, mu);
            
            % Compute Magnitude of velocity
            % Eps = V^2/2 - mu/r
            V = sqrt(2*(eps + mu./r));
            
            % Compute Tangential Velocity
            % h = r * vTangential
            vTangential = h./r;
            
            % Compute radial velocity sign
            vRsign = (-1).^ ( ta>pi );
            
            % Compute Radial Velocity
            vRadial = vRsign*sqrt(V.^2 - vTangential.^2);
            
        end
    end
end