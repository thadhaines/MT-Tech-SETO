% Example of miniwecc with AGC and selectable area config and VTS
clear all; close all; clc

%% Add pst path to MATLAB
% generate relative path generically
folderDepth = 3; % depth of current directory from main PST directory
pstVer =   'PSTv4';
pathParts = strsplit(pwd, filesep);
PSTpath = pathParts(1);

for pNdx = 2:max(size(pathParts))-folderDepth
    PSTpath = [char(PSTpath), filesep, char(pathParts(pNdx))];
end
PSTpath = [char(PSTpath), filesep, pstVer, filesep];

addpath(PSTpath)
save PSTpath.mat PSTpath pstVer
clear folderDepth pathParts pNdx PSTpath scenario

%% Prepare simulation files
clear all; close all; clc
load PSTpath.mat

areaSetting = 2;    % 1 for PSLTDSim Areas, 2 for EIA Areas
tripGen = 1;        % 1 for gen trip scenario, else load step
enableVTS = 1;      % 1 for VTS integration

% ensure file name accuracy
if enableVTS
    solnMethod = '_VTS';
else
    solnMethod = '_FTS';
end

% select which data file
if areaSetting
    copyfile('d_minniWECC_V3C_C3_6_C_NoFault_AGC_VTS_PSLTDSim.m',[PSTpath 'DataFile.m']); % PSLTDSim 3 Area
    areaConfig = 'PSLTDSim';
else
    copyfile('d_minniWECC_V3C_C3_6_C_NoFault_AGC_VTS_EIA.m',[PSTpath 'DataFile.m']); % 'EIA' 3 Area
    areaConfig = 'EIA';
end

% handle perturbance
if tripGen
    %Use gen trips
    if strcmp(pstVer, 'PSTv4')
        copyfile('mac_trip_logic_Gen_1_13_G.m', [PSTpath 'mac_trip_logic.m']);
    else
        error('Wrong PST Version')
    end
    pert = '_trip';
else
    % load step
    if strcmp(pstVer, 'PSTv4')
        copyfile('ml_sig_Step.m', [PSTpath 'ml_sig.m']);
    else
        error('Wrong PST Version')
    end
    pert = '_step';
end

%% Run Simulation
livePlotFlag = 0; % not ideal for extended term sim - may cause crashes/freezing

if strcmp(pstVer , 'PSTv4')
    s_simu
else
    s_simu_Batch %Run PST
end

%% Clean up modulation file alterations.
copyfile([PSTpath 'mac_trip_logic_ORIG.m'], [PSTpath 'mac_trip_logic.m']);
copyfile([PSTpath 'ml_sig_ORIG.m'],[PSTpath 'ml_sig.m']);

%% Save cleaned output data
fileSaveName = ['mw', areaConfig, pert, solnMethod,'.mat'];
save(fileSaveName); %Save simulation outputs

%% temp file clean up
delete('PSTpath.mat')

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

%% Plot tripped Gen derivatives
figure
n = find(g.sys.t<20);
plot(g.sys.t(n),g.mac.dmac_spd(1,n),'k')
hold on
plot(g.sys.t(n),g.mac.dmac_spd(7,n),'r')
plot(g.sys.t(n),g.mac.dmac_spd(13,n),'g','LineWidth',2)
grid on
ylabel('dGen Speed [PU]')
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
%xlim([0,20])
grid on
subplot(2,1,2)
plot(g.sys.t, g.sys.totH)
title('System Inertia')
%xlim([0,20])
ylabel('System Inertia [PU]')
xlabel('Time [seconds]')
grid on

%% time step plot
vts = zeros(g.vts.dataN-1,1);
for n=2:g.vts.dataN
    vts(n-1)= g.sys.t(n)-g.sys.t(n-1);
end

figure
stairs(g.sys.t(1:end-1),vts,'k','linewidth',1.25)
grid on
title({'Time Step Comparison'})
ylabel('Time Step Size [seconds]')
xlabel('Time [seconds]')

%% AGC capacity
figure
plot(g.sys.t, g.agc.agc(1).curGen/g.agc.agc(1).maxGen)
hold on
plot(g.sys.t, g.agc.agc(2).curGen/g.agc.agc(2).maxGen)
plot(g.sys.t, g.agc.agc(3).curGen/g.agc.agc(3).maxGen)
legend({'Area 1', 'Area 2', 'Area 3'}, 'location', 'best')
grid on
title({'Dispatched Area Generation'})
ylabel('Current Generation [PU]')
xlabel('Time [seconds]')
