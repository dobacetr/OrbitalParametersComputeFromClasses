classdef ClohessyWiltshire
    methods (Static)
        function Phirr = ComputePhirr(n, t)
            % Initialize
            Phirr = zeros(3, 'like',t);

            Phirr(1,1) = 4 - 3 * cos(n*t);

            Phirr(2,1) = 6*( sin(n*t)-n*t );
            Phirr(2,2) = 1;

            Phirr(3,3) = cos(n*t);
        end

        function Phirv = ComputePhirv(n, t)
            % Initialize
            Phirv = zeros(3, 'like',t);

            Phirv(1,1) = sin(n*t)/n;
            Phirv(1,2) = 2*( 1 - cos(n*t) )/n;

            Phirv(2,1) = 2 * ( cos(n*t) - 1 )/n;
            Phirv(2,2) = ( 4 * sin(n*t) - 3*n*t )/n;

            Phirv(3,3) = sin(n*t)/n;
        end

        function Phivr = ComputePhivr(n, t)
            % Initialize
            Phivr = zeros(3, 'like',t);

            Phivr(1,1) = 3 * n * sin(n*t);

            Phivr(2,1) = 6 * n * ( cos(n*t) - 1 );
            
            Phivr(3,3) = -n * sin(n*t);
        end

        function Phivv = ComputePhivv(n,t)
            % Initialize
            Phivv = zeros(3, 'like',t);

            Phivv(1,1) = cos(n*t);
            Phivv(1,2) = 2*sin(n*t);

            Phivv(2,1) = - 2* sin(n*t);
            Phivv(2,2) = 4*cos(n*t)-3;

            Phivv(3,3) = cos(n*t);
        end

    end
end