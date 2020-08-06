% Example of two machine Pmech ramp with govs
% Experimentation with pmech input mod

clear all; close all; clc

%% Add pst path to MATLAB
% generate relative path generically
folderDepth = 2; % depth of current directory from main PST directory
pstVer = 'pstSETO';
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
delete([PSTpath 'DataFile.m']); % ensure batch datafile is cleared
copyfile('d2a_ExtendedTermVTS.m',[PSTpath 'DataFile.m']); % copy system data file to batch run location

% Handle turbine modulation file placement
delete([PSTpath 'mpm_sig.m']); 
copyfile('mpm_sig_sinRampSETO.m',[PSTpath 'mpm_sig.m']); 

% use defined models
copyfile([PSTpath 'pss2.m'],[PSTpath 'pss.m']);
copyfile([PSTpath, 'mac_sub_NEW2.m'],[PSTpath ,'mac_sub.m']); 

s_simu_BatchVTS %Run PST

%% Simulation variable cleanup
% Clear any varables that contain only zeros
varNames = who()'; % all variable names in workspace
clearedVars = {}; % cell to hold names of deleted 'all zero' variables

for vName = varNames
    try
    zeroTest = eval(sprintf('all(%s(:)==0)', vName{1})); % check if all zeros
    if zeroTest
        eval(sprintf('clear %s',vName{1}) ); % clear variable
        clearedVars{end+1} = vName{1}; % add name to cell for reference
    end
    catch ME
        % gets called for structs... (global g)
        disp(ME.message)
        disp(vName)
    end

end
clear varNames vName zeroTest

%% Save cleaned output data
save('pmRampVTS.mat'); %Save simulation outputs

%% Clean up modulation file alterations.
% turbine governor pref mod
load PSTpath.mat
delete([PSTpath 'mpm_sig.m']); % remove simulation specific ml_sig file
copyfile([PSTpath 'mpm_sig_ORIG.m'],[PSTpath 'mpm_sig.m']); % Replace original file
copyfile([PSTpath 'mac_sub_ORIG.m'],[PSTpath 'mac_sub.m']);
copyfile([PSTpath 'pss3.m'],[PSTpath 'pss.m']);
%% temp file clean up
delete('PSTpath.mat')
