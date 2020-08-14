% Example of miniwecc in VTS
%
clear all; close all; clc

%% Add pst path to MATLAB
% generate relative path generically
folderDepth = 3; % depth of current directory from main PST directory
pstVer =   'PSTv4';%  'pstV3p1'; % 'pstV2P3'; %  
pathParts = strsplit(pwd, filesep);
PSTpath = pathParts(1);

for pNdx = 2:max(size(pathParts))-folderDepth
    PSTpath = [char(PSTpath), filesep, char(pathParts(pNdx))];
end
PSTpath = [char(PSTpath), filesep, pstVer, filesep];

addpath(PSTpath)
addpath([PSTpath, 'test', filesep]) % to handle new functionized code
save PSTpath.mat PSTpath pstVer
clear folderDepth pathParts p Ndx PSTpath

%% Run nonlinear simulation and store results
clear all; close all; clc
load PSTpath.mat
delete([PSTpath 'DataFile.m']); % ensure batch datafile is cleared
copyfile('d_minniWECC_VTSFault.m',[PSTpath 'DataFile.m']); 

copyfile([PSTpath 'pss2.m'],[PSTpath 'pss.m']); % specify pss model
copyfile([PSTpath 'mac_sub_NEW2.m'],[PSTpath 'mac_sub.m']); % specify machine model
copyfile([PSTpath 'livePlot_1.m'],[PSTpath 'livePlot.m']); % specify plot operation

livePlotFlag = 1;
% s_simu_Batch
s_simu %Run PST <- this is the main file to look at for simulation workings

%% Restore original models
copyfile([PSTpath 'pss3.m'],[PSTpath 'pss.m']); % specify pss model
copyfile([PSTpath 'mac_sub_ORIG.m'],[PSTpath 'mac_sub.m']); % specify pss model
copyfile([PSTpath 'livePlot_ORIG.m'],[PSTpath 'livePlot.m']); % specify pss model


%% Save cleaned output data
save('VTSFault0.mat'); %Save simulation outputs

%% temp file clean up
delete('PSTpath.mat')