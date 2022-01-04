classdef OrbitPeriod
    methods (Static)
        % Given Semi-Major axis(a), gravitational parameter(mu)
        function T = ComputeFrom_a_mu(a, mu)
            T = 2*pi ./ sqrt( mu ) .* a .^(3/2);
        end
    end
end