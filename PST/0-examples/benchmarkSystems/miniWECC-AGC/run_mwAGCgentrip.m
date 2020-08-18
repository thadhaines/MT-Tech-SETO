% Example of miniwecc with generator tripping

clear all; close all; clc

%% Add pst path to MATLAB
% generate relative path generically
folderDepth = 3; % depth of current directory from main PST directory
pstVer =   'PSTv4'; % 'pstV3p1'; %   'pstV2P3'; % 'pstSETO';%
pathParts = strsplit(pwd, filesep);
PSTpath = pathParts(1);

for pNdx = 2:max(size(pathParts))-folderDepth
    PSTpath = [char(PSTpath), filesep, char(pathParts(pNdx))];
end
PSTpath = [char(PSTpath), filesep, pstVer, filesep];

addpath(PSTpath)
save PSTpath.mat PSTpath pstVer
clear folderDepth pathParts pNdx PSTpath scenario

%% Run nonlinear simulation and store results
clear all; close all; clc
load PSTpath.mat

copyfile('d_minniWECC_V3C_C3_6_C_NoFault_AGC.m',[PSTpath 'DataFile.m']);
copyfile([PSTpath 'livePlot_2.m'],[PSTpath 'livePlot.m']);

%Use gen trips
if strcmp(pstVer, 'pstSETO') || strcmp(pstVer, 'PSTv4')
    copyfile('mac_trip_logic_Gen_1_13_G.m', [PSTpath 'mac_trip_logic.m']);
else
    copyfile('mac_trip_logic_Gen_1_13.m', [PSTpath 'mac_trip_logic.m']);
end

% % load step
% if strcmp(pstVer, 'pstSETO') || strcmp(pstVer, 'PSTv4')
%     copyfile('ml_sig_Step.m', [PSTpath 'ml_sig.m']);
% else
%     % not handled.
% end

livePlotFlag = 1;

if strcmp(pstVer , 'PSTv4')
    s_simu
else
    s_simu_Batch %Run PST 
end

%% Clean up modulation file alterations.
copyfile([PSTpath 'mac_trip_logic_ORIG.m'], [PSTpath 'mac_trip_logic.m']);
copyfile([PSTpath 'livePlot_ORIG.m'],[PSTpath 'livePlot.m']);
copyfile([PSTpath 'ml_sig_ORIG.m'],[PSTpath 'ml_sig.m']);

%% Save cleaned output data
save('tripTest.mat'); %Save simulation outputs

%% temp file clean up
delete('PSTpath.mat')


%% Plot 
figure
n = find(g.sys.t<20);
plot(g.sys.t(n),g.mac.mac_spd(1,n),'k')
hold on
plot(g.sys.t(n),g.mac.mac_spd(7,n),'r')
plot(g.sys.t(n),g.mac.mac_spd(13,n),'g','LineWidth',2)
ylabel('Gen Speed [PU]')
xlabel('Time [seconds]')
grid on
legend('Gen 1','Gen 7','Gen 13','Location','Best')

%% Plot tripped Gen speed
figure
n = find(g.sys.t<20);
plot(g.sys.t(n),g.mac.mac_spd(1,n),'k')
hold on
plot(g.sys.t(n),g.mac.mac_spd(7,n),'r')
plot(g.sys.t(n),g.mac.mac_spd(13,n),'g','LineWidth',2)
grid on
ylabel('Gen Speed [PU]')
xlabel('Time [seconds]')
legend('Gen 1','Gen 7','Gen 13','Location','Best')
title('Select Generator Speeds')

%% plot inertia
figure 
subplot(2,1,1)
plot(g.sys.t, g.sys.aveF)
ylabel('Ave. Sys. Freq. [PU]')
xlabel('Time [seconds]')
title('Average Weighted System Frequency')
xlim([0,20])
grid on
subplot(2,1,2)
plot(g.sys.t, g.sys.totH)
title('System Inertia')
xlim([0,20])
ylabel('System Inertia [PU]')
xlabel('Time [seconds]')
grid on