%% Run plots to compare version output
% requires comparePlots folder to be added to MATLAB path.
close all; clear; clc

dataName = 'SVCnonLIN';
a = ['pstSETO', dataName];
b = ['pstV2P3', dataName];
c = ['pstV3p1', dataName];

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


% difference in saturation of exc_dc12 between version 2 and 3