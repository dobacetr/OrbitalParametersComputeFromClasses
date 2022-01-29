function Cis = dcm_rsp2intertial(Ome, inc, w, ta)
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
end