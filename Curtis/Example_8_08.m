% ˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜
% Example_8_08
% ˜˜˜˜˜˜˜˜˜˜˜˜
%
% This program uses Algorithm 8.2 to solve Example 8.8.
%
% mu - gravitational parameter of the sun (km^3/s^2)
% deg - conversion factor between degrees and radians
% pi - 3.1415926...
%
% planet_id - planet identifier:
% 1 = Mercury
% 2 = Venus
% 3 = Earth
% 4 = Mars
% 5 = Jupiter
% 6 = Saturn
% 7 = Uranus
% 8 = Neptune
% 9 = Pluto
% planet_name - name of the planet
%
% year - range: 1901 - 2099
% month - range: 1 - 12
% month_name - name of the month
% day - range: 1 - 31
% hour - range: 0 - 23
% minute - range: 0 - 60
% second - range: 0 - 60
%
% depart - [planet_id, year, month, day, hour, minute,
% second] at departure
% arrive - [planet_id, year, month, day, hour, minute,
% second] at arrival
%
% planet1 - [Rp1, Vp1, jd1]
% planet2 - [Rp2, Vp2, jd2]
% trajectory - [V1, V2]
%
% coe - orbital elements [h e RA incl w TA]
% where
% h = angular momentum (km^2/s)
% e = eccentricity
% RA = right ascension of the ascending
% node (rad)
% incl = inclination of the orbit (rad)
% w = argument of perigee (rad)
% TA = true anomaly (rad)
% a = semimajor axis (km)
%
% jd1, jd2 - Julian day numbers at departure and arrival
% tof - time of flight from planet 1 to planet 2
% (days)
%
% Rp1, Vp1 - state vector of planet 1 at departure
% (km, km/s)
% Rp2, Vp2 - state vector of planet 2 at arrival
% (km, km/s)
% R1, V1 - heliocentric state vector of spacecraft at
% departure (km, km/s)
% R2, V2 - heliocentric state vector of spacecraft at
% arrival (km, km/s)
%
% vinf1, vinf2 - hyperbolic excess velocities at departure
% and arrival (km/s)
%
% User M-functions required: interplanetary, coe_from_sv,
% month_planet_names
% ------------------------------------------------------------
clear
global mu
mu = 1.327124e11;
deg = pi/180;
%...Data for planet 1:
planet_id = 3; % (earth)
year = 1996;
month = 11;
day = 7;
hour = 0;
minute = 0;
second = 0;
%...
depart = [planet_id year month day hour minute second];
%...Data for planet 2:
planet_id = 4; % (Mars)
year = 1997;
month = 9;
day = 12;
hour = 0;
minute = 0;
second = 0;
%...
arrive = [planet_id year month day hour minute second];
[planet1, planet2, trajectory] = interplanetary ...
(depart, arrive);
R1 = planet1(1,1:3);
Vp1 = planet1(1,4:6);
jd1 = planet1(1,7);
R2 = planet2(1,1:3);
Vp2 = planet2(1,4:6);
jd2 = planet2(1,7);
V1 = trajectory(1,1:3);
V2 = trajectory(1,4:6);
tof = jd2 - jd1;
%...Use Algorithm 4.1 to find the orbital elements of the
% spacecraft trajectory based on [Rp1, V1]...
coe = coe_from_sv(R1, V1);
% ... and [R2, V2]
coe2 = coe_from_sv(R2, V2);
%...Equations 8.102 and 8.103:
vinf1 = V1 - Vp1;
vinf2 = V2 - Vp2;
%...Echo the input data and output the solution to
% the command window:
fprintf('---------------------------------------------------')
fprintf('\n Example 8.8')
fprintf('\n\n Departure:\n');
[month_name, planet_name] = month_planet_names(depart(3), ...
depart(1));
fprintf('\n Planet: %s', planet_name)
fprintf('\n Year : %g', depart(2))
fprintf('\n Month : %s', month_name)
fprintf('\n Day : %g', depart(4))
fprintf('\n Hour : %g', depart(5))
fprintf('\n Minute: %g', depart(6))
fprintf('\n Second: %g', depart(7))
fprintf('\n\n Julian day: %11.3f\n', jd1)
fprintf('\n Planet position vector (km) = [%g %g %g]', ...
R1(1), R1(2), R1(3))
fprintf('\n Magnitude = %g\n', norm(R1))
fprintf('\n Planet velocity (km/s) = [%g %g %g]', ...
Vp1(1), Vp1(2), Vp1(3))
fprintf('\n Magnitude = %g\n', norm(Vp1))
fprintf('\n Spacecraft velocity (km/s) = [%g %g %g]', ...
V1(1), V1(2), V1(3))
fprintf('\n Magnitude = %g\n', norm(V1))
fprintf('\n v-infinity at departure (km/s) = [%g %g %g]', ...
vinf1(1), vinf1(2), vinf1(3))
fprintf('\n Magnitude = %g\n', norm(vinf1))
fprintf('\n\n Time of flight = %g days\n', tof)
fprintf('\n\n Arrival:\n');
[month_name, planet_name] = month_planet_names(arrive(3), ...
arrive(1));
fprintf('\n Planet: %s', planet_name)
fprintf('\n Year : %g', arrive(2))
fprintf('\n Month : %s', month_name)
fprintf('\n Day : %g', arrive(4))
fprintf('\n Hour : %g', arrive(5))
fprintf('\n Minute: %g', arrive(6))
fprintf('\n Second: %g', arrive(7))
fprintf('\n\n Julian day: %11.3f\n', jd2)
fprintf('\n Planet position vector (km) = [%g %g %g]', ...
R2(1), R2(2), R2(3))
fprintf('\n Magnitude = %g\n', norm(R1))
fprintf('\n Planet velocity (km/s) = [%g %g %g]', ...
Vp2(1), Vp2(2), Vp2(3))
fprintf('\n Magnitude = %g\n', norm(Vp2))
fprintf('\n Spacecraft Velocity (km/s) = [%g %g %g]', ...
V2(1), V2(2), V2(3))
fprintf('\n Magnitude = %g\n', norm(V2))
fprintf('\n v-infinity at arrival (km/s) = [%g %g %g]', ...
vinf2(1), vinf2(2), vinf2(3))
fprintf('\n Magnitude = %g', norm(vinf2))
fprintf('\n\n\n Orbital elements of flight trajectory:\n')
fprintf('\n Angular momentum (km^2/s) = %g', coe(1))
fprintf('\n Eccentricity = %g', coe(2))
fprintf('\n Right ascension of the ascending node')
fprintf(' (deg) = %g', coe(3)/deg)
fprintf('\n Inclination to the ecliptic (deg) = %g', ...
coe(4)/deg)
fprintf('\n Argument of perihelion (deg) = %g', ...
coe(5)/deg)
fprintf('\n True anomaly at departure (deg) = %g', ...
coe(6)/deg)
fprintf('\n True anomaly at arrival (deg) = %g\n', ...
coe2(6)/deg)
fprintf('\n Semimajor axis (km) = %g', coe(7))
if coe(2) < 1
fprintf('\n Period (days) = %g', ...
2*pi/sqrt(mu)*coe(7)^1.5/24/3600)
end
fprintf('\n-----------------------------------------------\n')
% ˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜
% Output from Example_8_08
% -----------------------------------------------------
% Example 8.8
% Departure:
% Planet: Earth
% Year : 1996
% Month : November
% Day : 7
% Hour : 0
% Minute: 0
% Second: 0
% Julian day: 2450394.500
% Planet position vector (km) = [1.04994e+08 1.04655e+08 988.331]
% Magnitude = 1.48244e+08
% Planet velocity (km/s) = [-21.515 20.9865 0.000132284]
% Magnitude = 30.0554
% Spacecraft velocity (km/s) = [-24.4282 21.7819 0.948049]
% Magnitude = 32.7427
% v-infinity at departure (km/s) = [-2.91321 0.79542 0.947917]
% Magnitude = 3.16513
% Time of flight = 309 days
% Arrival:
% Planet: Mars
% Year : 1997
% Month : September
% Day : 12
% Hour : 0
% Minute: 0
% Second: 0
% Julian day: 2450703.500
% Planet position vector (km) = [-2.08329e+07 -2.18404e+08 -4.06287e+06]
% Magnitude = 1.48244e+08
% Planet velocity (km/s) = [25.0386 -0.220288 -0.620623]
% Magnitude = 25.0472
% D.18 Algorithm 8.2: calculation of the spacecraft trajectory 655
% Spacecraft Velocity (km/s) = [22.1581 -0.19666 -0.457847]
% Magnitude = 22.1637
% v-infinity at arrival (km/s) = [-2.88049 0.023628 0.162776]
% Magnitude = 2.88518
% Orbital elements of flight trajectory:
% Angular momentum (km^2/s) = 4.84554e+09
% Eccentricity = 0.205785
% Right ascension of the ascending node (deg) = 44.8942
% Inclination to the ecliptic (deg) = 1.6621
% Argument of perihelion (deg) = 19.9738
% True anomaly at departure (deg) = 340.039
% True anomaly at arrival (deg) = 199.695
% Semimajor axis (km) = 1.84742e+08
% Period (days) = 501.254
% -------------------------------------------------------------