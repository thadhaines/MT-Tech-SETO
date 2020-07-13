%% Run plots to compare version output
% requires comparePlots folder to be added to MATLAB path.
close all; clear; clc

dataName = 'INDnonLIN';
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

% %% Compare svc_sig, B_cv...
% compareB_cv( a, c, printFigs )
% compareB_cv( a, b, printFigs )
% compareB_cv( c, b, printFigs )


% difference in angle due to mac_ind model - two versions exist