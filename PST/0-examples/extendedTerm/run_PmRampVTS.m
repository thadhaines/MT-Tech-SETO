% Example of two machine Pmech ramp with govs
% Experimentation with pmech input mod

clear all; close all; clc

%% Add pst path to MATLAB
% generate relative path generically
folderDepth = 2; % depth of current directory from main PST directory
pstVer = 'PSTv4';
pathParts = strsplit(pwd, filesep);
PSTpath = pathParts(1);

for pNdx = 2:max(size(pathParts))-folderDepth
    PSTpath = [char(PSTpath), filesep, char(pathParts(pNdx))];
end
PSTpath = [char(PSTpath), filesep, pstVer, filesep];

addpath(PSTpath)

save PSTpath.mat PSTpath
clear folderDepth pathParts pNdx PSTpath

%% Run nonlinear simulation and store results
clear all; close all; clc
load PSTpath.mat

copyfile('d2a_ExtendedTermVTS.m',[PSTpath 'DataFile.m']); % copy system data file to batch run location

% Handle turbine modulation file placement
copyfile('mpm_sig_sinRampSETO.m',[PSTpath 'mpm_sig.m']); 

% use defined models
copyfile([PSTpath 'pss3.m'],[PSTpath 'pss.m']);
copyfile([PSTpath, 'mac_sub_NEW2.m'],[PSTpath ,'mac_sub.m']); 

s_simu %Run PST

%% Save cleaned output data
save('pmRampVTS.mat'); %Save simulation outputs

%% Clean up modulation file alterations.
% turbine governor pref mod
load PSTpath.mat

copyfile([PSTpath 'pss2.m'],[PSTpath 'pss.m']);
copyfile([PSTpath 'mpm_sig_ORIG.m'],[PSTpath 'mpm_sig.m']); % Replace original file
copyfile([PSTpath 'mac_sub_ORIG.m'],[PSTpath 'mac_sub.m']);

%% temp file clean up
delete('PSTpath.mat')
