% This class contains methods for computing OrbitRadius
classdef OrbitPosition
    
    
    %% Generic Functions
    methods (Static)
        
        function rVec = ComputeFrom_h_e_inc_Ome_w_ta_mu(h, e, inc, Ome, w, ta, mu)
            r = OrbitRadius.ComputeFrom_h_e_ta_mu(h, e, ta, mu);
            
            rVec = OrbitPosition.ComputeFrom_inc_Ome_w_ta_r(inc, Ome, w, ta, r);
        end
        
        % Given orbit inc, Ome, w, ta, r; Compute 3D position
        function rVec = ComputeFrom_inc_Ome_w_ta_r(inc, Ome, w, ta, r)
            % Compute Cip, DCM Perifocal To Inertial 
            Cip = permute(dcm_eci2perifocal(Ome, inc, w), [2,1,3]);
            
            % Compute DCM rsp(radial, tangential , perpendicular) to perifocal
            cta = cos(ta);
            sta = sin(ta);
            szta = length(ta);
            Cpr = reshape(...
                [...
                cta;
                sta;
                zeros(1, szta);
                -sta;
                cta;
                zeros(3, szta);
                ones(1, szta);
                ], ...
                [3, 3, szta]);
            
            % Compute DCM rsp to inertial
            Cis = Mmult( Cip, Cpr );
            
            % Position in srp is [r;0;0]
            rVec = MVmult(Cis, [r; zeros(2, length(r))]);
            
        end
    end
end