function JD = ComputeJDFromEpoch(epochYear, epochDay, epochMonth)
% epoch year is last 2 digits
% it is considered between 1957 and 2056 for now
assert(epochYear<100, 'epochYear must be the last 2 digits of year.');
if(epochYear<57)
    epochYear = epochYear + 100;
end
epochYear = epochYear + 1900;

if nargin < 3
    monthLengths = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31];

    % Compute leap day
    % every 4th year, but not 100th, but 400th is a year with leap day
    if mod(epochYear,4)==0 && ( mod(epochYear,100)~=0 || mod(epochYear,400)==0)
        monthLengths(2) = monthLengths(2) + 1;
    end
    monthBreaks = cumsum(monthLengths);

    % Search where the day is at
    epochMonth = 1;
    while monthBreaks(epochMonth) < epochDay && epochMonth < 12
        epochMonth = epochMonth + 1;
    end
    
    % subtract beginning of the month
    if epochMonth>1
        epochDay = epochDay - monthBreaks(epochMonth-1);
    end
end

% Compute hour from day's fractional component and dispose the fractional
epochHour = (epochDay - floor(epochDay))*24;
epochDay  = floor(epochDay);

% Compute minute from hour's fractional component and dispose the fractional
epochMinute = ( epochHour - floor(epochHour) ) * 60;
epochHour = floor(epochHour);

% Compute second from minute's fractional component and dispose the fractional
epochSecond = ( epochMinute - floor(epochMinute) ) * 60;
epochMinute = floor(epochMinute);

JD = ComputeJD(epochYear,epochMonth,epochDay,epochHour,epochMinute,epochSecond);
end