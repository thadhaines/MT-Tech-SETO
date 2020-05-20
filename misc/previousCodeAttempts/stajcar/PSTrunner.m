% Matt Stajcar
% PST runner

clear
close all
clc

addpath('D:\Users\mstajcar\Documents\2014-2015\Summer\Simulation Project\pstV2')

s_simu

%% plot
figure
plot(t,mac_spd(1,:)*60,'b',t,mac_spd(2,:)*60,'r')
grid on

figure
plot(t,mac_ang(1,:)*180/pi,'b',t,mac_ang(2,:)*180/pi,'r')
grid on

tPST = t;
w1PST = mac_spd(1,:);
w2PST = mac_spd(2,:);

ang1PST = mac_ang(1,:);
ang2PST = mac_ang(2,:);

save PST tPST w1PST w2PST ang1PST ang2PST

