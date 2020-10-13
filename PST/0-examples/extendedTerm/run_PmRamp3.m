% Example of two machine Pmech ramp with govs
% Experimentation with pmech input mod

clear all; close all; clc

%% Add pst path to MATLAB
% generate relative path generically
folderDepth = 2; % depth of current directory from main PST directory
pstVer = 'pstV3p1';

pathParts = strsplit(pwd, filesep);
PSTpath = pathParts(1);

for pNdx = 2:max(size(pathParts))-folderDepth
    PSTpath = [char(PSTpath), filesep, char(pathParts(pNdx))];
end
PSTpath = [char(PSTpath), filesep, pstVer, filesep];

addpath(PSTpath)
addpath([PSTpath,'test',filesep]) % for test VTS functions
save PSTpath.mat PSTpath
clear folderDepth pathParts pNdx PSTpath

%% Run nonlinear simulation and store results
clear all; close all; clc
load PSTpath.mat

copyfile('d2a_ExtendedTerm.m',[PSTpath 'DataFile.m']); % copy system data file to batch run location

% Handle turbine modulation file placement etc...
copyfile('mpm_sig_sinRamp3.m',[PSTpath 'mpm_sig.m']); % copy simulation specific data file to batch run location

s_simu_Batch %Run PST

%% Save cleaned output data
save('pmRamp3FTS.mat'); %Save simulation outputs

%% Clean up modulation file alterations.
% turbine governor pref mod
%load PSTpath.mat
copyfile([PSTpath 'mpm_sig_ORIG.m'],[PSTpath 'mpm_sig.m']); % Replace original file

%% temp file clean up
delete('PSTpath.mat')
delete('sim_fle.mat')
