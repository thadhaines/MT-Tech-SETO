clear all; close all; clc
PSTpath = 'C:\Users\dtrudnowski\Documents\PST\pstV2p3\';
addpath(PSTpath)
delete([PSTpath 'DataFile.m']); copyfile('d_minniWECC_V3C_C3_6_C_NoFault.m',[PSTpath 'DataFile.m']); %System data file
delete([PSTpath 'mac_trip_logic.m']); copyfile('mac_trip_logic_Gen_1_13.m',[PSTpath 'mac_trip_logic.m']); %tripping files
s_simu_Batch %Run PST
delete([PSTpath 'mac_trip_logic.m']); copyfile([PSTpath 'mac_trip_logic_ORIG.m'],[PSTpath 'mac_trip_logic.m']); %reset tripping file

%% Plot 
figure
n = find(t<20);
plot(t(n),mac_spd(1,n),'k',t(n),mac_spd(7,n),'r',t(n),mac_spd(13,n),'g','LineWidth',2)
yabel('gen spd (pu)')
xlabel('Time (s)')
legend('Gen 1','Gen 7','Gen 13','Location','Best')