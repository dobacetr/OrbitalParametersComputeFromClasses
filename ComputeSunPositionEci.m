function rSunEci = ComputeSunPositionEci(JD)
    % Taken from
    % https://astronomy.stackexchange.com/questions/28802/calculating-the-sun-s-position-in-eci
    % Tested against Curtis' code for Earth's Position around Sun Orbit.
    % 0.2 degrees different unit vectors were observed. This might be more
    % accurate as the formet was designed for interplanetary problems.
    % J000 = ComputeJDFromEpoch(0, 1.5);
    J2000 = 2451545; % Jan 1, 2000, 12:00
    
    d = JD - J2000;

    % Calculate parameters
    L = 280.4606184 + (36000.77005361 / 36525) * d;% a.k.a. mean longitude, in degrees
    g = 357.5277233 + (35999.05034 / 36525 * d);% a.k.a. mean anomaly, in degrees
    p = L + 1.914666471 * sin(g * pi / 180) + 0.918994643 * sin(2*g * pi / 180);% a.k.a. ecliptic longitude lambda, in degrees
    q = 23.43929 - ((46.8093/3600) * (d / 36525)) - 0.0001831/3600*(d/36525)^2;% a.k.a. obliquity of ecliptic plane epsilon, in degrees

    % 2. Calculate unit directional vector in ECI coordinates
    u_x = cos(p * pi / 180);%
    u_y = cos(q * pi / 180) * sin(p * pi / 180);%
    u_v = sin(q * pi / 180) * sin(p * pi / 180);%
    u = [u_x;u_y;u_v];

    % 3. Calculate distance to sun and scale the unit vector
    a = 1.000140612 - 0.016708617 * cos(g * pi / 180) - 0.000139589 * cos(2 * g * pi / 180);% a.k.a distance from Earth's center to Sun's center in astronomical units (AU)
    m = a * 149597870700;% a.k.a. center-to-center distance from Earth to Sun in meters
    rSunEci = m * u;% (or a * u_v) is distance to sun in meters (or in AU)
end