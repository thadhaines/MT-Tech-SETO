%% Run plots to compare version output
% requires comparePlots folder to be added to MATLAB path.
close all; clear; clc

a = 'pstSETOd19new';
b = 'pstV2P3d19new';
c = 'pstV3p1d19new';

printFigs = 0;
%% Bus Voltage
compareBus_V( a, c, printFigs )
compareBus_V( a, b, printFigs )
compareBus_V( c, b, printFigs )

%% Bus angle
compareBus_Angle( a, c, printFigs )
compareBus_Angle( a, b, printFigs )
compareBus_Angle( c, b, printFigs )

%% Machine Speed
compareMac_Spd( a, c, printFigs )
compareMac_Spd( a, b, printFigs )
compareMac_Spd( c, b, printFigs )

%differences in machine speed leads to voltage differences