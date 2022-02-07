classdef OrbitTime
    methods (Static)
        function t = ComputeTimeFromPerigee_Elliptic(T,e,ta)
            atanTerm = atan( sqrt( (1-e)/(1+e) )*tan(ta/2) );
            t = T/(2*pi)*...
                (...
                2*atanTerm ...
                - e*sin( 2*atanTerm ) ...
                );
        end
        
        function t = ComputeTimeFromPerigee_Hyperbolic(h,mu,e,ta)
            atanhTerm = atanh( sqrt( (e-1)/(1+e) )*tan(ta/2) );
            
            t = h^3/( mu^2*( e^2-1 )^(3/2) )*...
                ( e*sinh(2*atanhTerm) -2*atanhTerm );
        end
    end
end