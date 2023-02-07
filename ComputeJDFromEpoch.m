function JD = ComputeJDFromEpoch(epochYear, epochDay, epochMonth)
% J2000 : 1 Jan 2000, 12:00 is (0, 1.5) or (0, 1.5, 0)
% 7.2.2023 16:00 is (23, 7+16/24, 2)


% epoch year is last 2 digits
% it is considered between 1957 and 2056 for now
assert(epochYear<100, 'epochYear must be the last 2 digits of year.');
if(epochYear<57)
    epochYear = epochYear + 100;
end
epochYear = epochYear + 1900;

if nargin < 3
    [epochMonth, epochDay, epochHour, epochMinute, epochSecond] = EpochDay2EpochMonthDayHourMinuteSecond(epochYear, epochDay);
else
    [epochMonth, epochDay, epochHour, epochMinute, epochSecond] = EpochDay2EpochMonthDayHourMinuteSecond(epochYear, epochDay, epochMonth);
end

JD = ComputeJD(epochYear,epochMonth,epochDay,epochHour,epochMinute,epochSecond);
end