function rSunEci = ComputeSunPositionEci(JD)
    % Taken from
    % https://astronomy.stackexchange.com/questions/28802/calculating-the-sun-s-position-in-eci
    % Tested against Curtis' code for Earth's Position around Sun Orbit.
    % 0.2 degrees different unit vectors were observed. This might be more
    % accurate as the former was designed for interplanetary problems.
    % J000 = ComputeJDFromEpoch(0, 1.5);
    J2000 = 2451545; % Jan 1, 2000, 12:00
    
    if size(JD, 1) > 1
        JD = reshape(JD, 1, []);
    end
    
    d = JD - J2000;

    % Calculate parameters
    L = (pi/180)*280.4606184 + ((pi/180)*36000.77005361 / 36525) * d;% a.k.a. mean longitude, in rad
    g = (pi/180)*357.5277233 + ((pi/180)*35999.05034 / 36525 )* d;% a.k.a. mean anomaly, in rad
    p = L + (pi/180)*1.914666471 * sin(g) + (pi/180)*0.918994643 * sin(2*g);% a.k.a. ecliptic longitude lambda, in rad
    q = (pi/180)*23.43929 - ((pi/180)*46.8093/(3600*36525)) * d - ((pi/180)*0.0001831/(3600*36525^2))*(d).^2;% a.k.a. obliquity of ecliptic plane epsilon, in rad

    % 2. Calculate unit directional vector in ECI coordinates
    sP = sin(p);
    u = [...
        cos(p);
        cos(q) .* sP;
        sin(q) .* sP;
        ];

    % 3. Calculate distance to sun and scale the unit vector
    a = 1.000140612 - 0.016708617 * cos(g) - 0.000139589 * cos(2 * g);% a.k.a distance from Earth's center to Sun's center in astronomical units (AU)
    m = a * 149597870700;% a.k.a. center-to-center distance from Earth to Sun in meters
    if coder.target('MATLAB')
        rSunEci = m .* u;% (or a * u_v) is distance to sun in meters (or in AU)
    else
        rSunEci = repmat(m, 3, 1) .* u;% (or a * u_v) is distance to sun in meters (or in AU)
    end
end