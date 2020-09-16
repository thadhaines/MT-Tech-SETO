% Example of two machin load step with govs
% smaller system used to more easily study inner workings of PST 

% Experimentation with load modulation and gov
% Load step on bus 2 of +0.25 PU at t=1 where lmod TS = 0.1 second

clear all; close all; clc

%% Add pst path to MATLAB
% generate relative path generically
folderDepth = 3; % depth of current directory from main PST directory
pstVer = 'PSTv4'; % 'pstSETO';  'pstV2P3'; % 'pstV3p1'; %
pathParts = strsplit(pwd, filesep);
PSTpath = pathParts(1);

for pNdx = 2:max(size(pathParts))-folderDepth
    PSTpath = [char(PSTpath), filesep, char(pathParts(pNdx))];
end
PSTpath = [char(PSTpath), filesep, pstVer, filesep];

addpath(PSTpath)
save PSTpath.mat PSTpath pstVer
clear folderDepth pathParts pNdx PSTpath

%% Run nonlinear simulation and store results
clear all; close all; clc
load PSTpath.mat

copyfile('d_2machineNOGov.m',[PSTpath 'DataFile.m']); % copy system data file to batch run location

% copy simulation specific data file to batch run location
if strcmp('PSTv4', pstVer) || strcmp('pstSETO', pstVer)
    copyfile('ml_sig_loadStepG.m',[PSTpath 'ml_sig.m']);
else
    copyfile('ml_sig_loadStep.m',[PSTpath 'ml_sig.m']);
end


if strcmp('PSTv4', pstVer)
    s_simu
else
    s_simu_Batch %Run PST
end

%% Clean up modulation file alterations.
copyfile([PSTpath 'ml_sig_ORIG.m'],[PSTpath 'ml_sig.m']);

%% Save cleaned output data
save('twoMachineNOGov.mat'); %Save simulation outputs

%% temp file clean up
delete('PSTpath.mat')
%%
figure
if strcmp('PSTv4', pstVer) || strcmp('pstSETO', pstVer)
    plot(g.sys.t, g.mac.mac_spd)
else
    plot(t, mac_spd)
end

xlabel('Time [seconds]')
ylabel('Machine Speed [PU]')
title('Machine Speeds')
grid on
