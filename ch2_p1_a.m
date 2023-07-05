% Volume =  pi*( ((30/2)^2*25) + ((65/2)^2*40) + ((70/2)^2*175) + ((90/2)^2*55) + ((99.6/2)^2*120) + ((90/2)^2*625) + ((70/2)^2*35) + ((50/2)^2*23) + ((40/2)^2*52))
% Len = 1150 
% syms d
% Diameter = double(abs(solve(Volume==pi*(d/2)^2*Len,d)))

% Volume =  pi*( ((30/2)^2*25) + ((65/2)^2*(32.5-25)) )
% Len = 25+(32.5-25)
% syms d
% Diameter = double(abs(solve(Volume==pi*(d/2)^2*Len,d)))

% Volume =  pi*( ((70/2)^2*35) + ((50/2)^2*23) + ((40/2)^2*52) )
% Len = 110
% syms d
% Diameter = double(abs(solve(Volume==pi*(d/2)^2*Len,d)))

% Volume =  pi*( ((70/2)^2*35) + ((50/2)^2*(23/2)) )
% Len = 46.5
% syms d
% Diameter = double(abs(solve(Volume==pi*(d/2)^2*Len,d)))

% Volume =  pi*( ((65.62/2)^2*46.5) + ((90/2)^2*(31.25)) )
% Len = 46.5+31.25
% syms d
% Diameter = double(abs(solve(Volume==pi*(d/2)^2*Len,d)))

% Volume =  pi*( ((70/2)^2*35) + ((90/2)^2*20) + ((110/2)^2*25) )
% Len = 80
% syms d
% Diameter = double(abs(solve(Volume==pi*(d/2)^2*Len,d)))

% Volume =  pi*( ((50/2)^2*11.5) + ((40/2)^2*(31.75-11.5)) )
% Len = 31.75
% syms d
% Diameter = double(abs(solve(Volume==pi*(d/2)^2*Len,d)))

% VolumeTondo = pi*(65/2)^2*40
% VolumeTrench = 2*(9*2.5*32.5)
% Volume = VolumeTondo - VolumeTrench
% Len = 40
% syms d
% Diameter = double(abs(solve(Volume==pi*(d/2)^2*Len,d)))

VolumeTondo = pi*(40/2)^2*52
VolumeTrench = 12*5*46
Volume = VolumeTondo - VolumeTrench
Len = 52
syms d
Diameter = double(abs(solve(Volume==pi*(d/2)^2*Len,d)))





