function RV = laplace_orbit_fit(lat, lon, alt, T, AZI_ELE)
    constants;
    for (i = 1:numel(T))
        dv = datevec(T(i));
        [RA_DEC(1,i),RA_DEC(2,i)] = horizontal_to_equatorial(AZI_ELE(1,i), ...
                AZI_ELE(2,i), lat, lon, dv(1), dv(2), dv(3), dv(4), dv(5), dv(6));
    end
    L = [cos(RA_DEC(1,:)).*cos(RA_DEC(2,:));
         sin(RA_DEC(1,:)).*cos(RA_DEC(2,:)); ...
         sin(RA_DEC(2,:))];
    LDOT1 = ((T(2)-T(3)) / ((T(1)-T(2))*(T(1)-T(3))) * L(:,1) + ...
            (2*T(2)-T(1)-T(3)) / ((T(2)-T(1))*(T(2)-T(3))) * L(:,2) + ...
            (T(2)-T(1)) / ((T(3)-T(1))*(T(3)-T(2))) * L(:,3))/86400;
    LDOT2 =(2 / ((T(1)-T(2))*(T(1)-T(3))) * L(:,1) + ...
            2 / ((T(2)-T(1))*(T(2)-T(3))) * L(:,2) + ...
            2 / ((T(3)-T(1))*(T(3)-T(2))) * L(:,3))/(86400*86400);
    dv = datevec(T(2));
    [RVA,~] = geodetic_to_ECI(lat, lon, alt, dv(1), dv(2), dv(3), dv(4), dv(5), dv(6));
    RVA = RVA/SMA_E;
    R = norm(RVA(1:3));
    N = dot(L(:,2), RVA(1:3));
    D  = det(horzcat(L(:,2), LDOT1, LDOT2))*2;
    D1 = det(horzcat(L(:,2), LDOT1, RVA(7:9)));
    D2 = det(horzcat(L(:,2), LDOT1, RVA(1:3)));
    D3 = det(horzcat(L(:,2), RVA(7:9), LDOT2));
    D4 = det(horzcat(L(:,2), RVA(1:3), LDOT2));
    c8 = D*D;
    c6 = ((4*N*D - 4*D1)*D1 - D*D*R*R);
    c3 = (4*MU_E_ER*D2*(N*D - 2*D1));
    c0 = -4*MU_E_ER*MU_E_ER*D2*D2;
    Z = roots([c8 0 c6 0 0 c3 0 0 c0]);
    for (i = 1:numel(Z))
        if (isreal(Z(i)) && Z(i) > 0)
            r = Z(i);
            break;
        end
    end
    rrr = r*r*r;
    rho = -2*(D1/D) - 2*(MU_E_ER/rrr)*(D2/D);
    rhodot = -(D3/D) - (MU_E_ER/rrr)*(D4/D);
    RV = vertcat(rho*L(:,2) + RVA(1:3), ...
            rhodot*L(:,2) + rho*LDOT1 + RVA(4:6))*SMA_E;
end
