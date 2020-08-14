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
delete([PSTpath 'DataFile.m']); % ensure batch datafile is cleared

copyfile('d_minniWECC_V3C_C3_6_C_NoFault.m',[PSTpath 'DataFile.m']);
copyfile([PSTpath 'livePlot_1.m'],[PSTpath 'livePlot.m']);

if strcmp(pstVer, 'pstSETO') || strcmp(pstVer, 'PSTv4') 
copyfile('mac_trip_logic_Gen_1_13_G.m', [PSTpath 'mac_trip_logic.m']);
else
copyfile('mac_trip_logic_Gen_1_13.m', [PSTpath 'mac_trip_logic.m']);
end

livePlotFlag = 1;

if strcmp(pstVer , 'PSTv4')
    s_simu
else
s_simu_Batch %Run PST <- this is the main file to look at for simulation workings
end

%% Clean up modulation file alterations.
copyfile([PSTpath 'mac_trip_logic_ORIG.m'], [PSTpath 'mac_trip_logic.m']);
copyfile([PSTpath 'livePlot_ORIG.m'],[PSTpath 'livePlot.m']);

%% Save cleaned output data
save('test.mat'); %Save simulation outputs

%% temp file clean up
delete('PSTpath.mat')