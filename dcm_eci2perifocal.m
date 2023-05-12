function C_eci2perifocal = dcm_eci2perifocal(OmeRad, incRad, argPeriRad)

    % Inputs:
    % OmeRad is right ascension of the ascending node in radians.
    % incRad is the inclination in radians.
    % argPeriRad is the argument of periapsis in radians.

    % From ECI to perifocal frame
    % Turn by RightAscensionOfTheAscendingNode around z
    %   Now, x looks at the ascending node
    % Turn by Inclination around x
    %   Now, we x and y are in the orbit plane
    % Turn by ArgumentOfPerigee around z
    %   Now, x looks at the Perigee
    
    cOme = cos(OmeRad);
    sOme = sin(OmeRad);
    
    cInc = cos(incRad);
    sInc = sin(incRad);
    
    cArg = cos(argPeriRad);
    sArg = sin(argPeriRad);
    
    szOme = length(OmeRad);
    szInc = length(incRad);
    szArg = length(argPeriRad);
    
%     % Craft C1 such that it is 3x3xN, C1 rotates coordinate system around
%     % positive z
%     C1 = reshape([...
%         cOme;
%         -sOme;
%         zeros(1, szOme);
%         sOme;
%         cOme;
%         zeros(3, szOme);
%         ones(1, szOme);], 3, 3, szOme);
%     
%     % Craft C2 such that it is 3x3xN, C2 rotates coordinate system around
%     % positive x
%     C2 = reshape([...
%         ones(1, szInc);
%         zeros(3, szInc);
%         cInc;
%         -sInc;
%         zeros(1, szInc);
%         sInc;
%         cInc], 3, 3, szInc);
%     
%     % Craft C3 such that it is 3x3xN, C3 rotates coordinate system around
%     % positive z
%     C3 = reshape([...
%         cArg;
%         -sArg;
%         zeros(1, szArg);
%         sArg;
%         cArg;
%         zeros(3, szArg);
%         ones(1, szArg);], 3, 3, szArg);
%     
%     % Compute the final
%     C_eci2perifocal = Mmult(Mmult(C3, C2) , C1);
%     C_eci2perifocal = DCM_angle(3, ome) * DCM_angle(1, inc) * DCM_angle(3, Ome)

 
C_eci2perifocal = reshape([...
    cArg.*cOme - cInc.*sArg.*sOme;
    - cOme.*sArg - cArg.*cInc.*sOme;
    sInc.*sOme;
    cArg.*sOme + cInc.*cOme.*sArg;
    cArg.*cInc.*cOme - sArg.*sOme;
    -cOme.*sInc;
    sArg.*sInc;
    cArg.*sInc;
    cInc], 3, 3, max([szOme, szInc, szArg]));
 
