function JD = ComputeJD(year,month,day,hour,minute,second)

% Ignore leap seconds
JD = 1721013.5 + 367*year - fix( 7/4*(year + fix((month+9)/12)) ) + fix(275*month/9) + day + (60*hour + minute + second/60)/1440;
% % Another Formula from the internet
% if month < 3
%     year = year-1;
%     month = month + 12;
% end
% 
% A = floor(year/100);
% B = floor(A/4);
% C = 2 - A + B;
% E = floor(365.25*(year+4716));
% F = floor(30.6001*(month+1));
% JD = C+day+E+F-1524.5;
% % add the hour minute second to fix the missing components
% JD = JD + (hour+(minute+second/60)/60)/24;odt