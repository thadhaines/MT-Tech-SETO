%% File to New England 39 bus system
clear all; close all; clc
caseName = 'datanp48';

%% Add pst path to MATLAB
% generate relative path generically
folderDepth = 3; % depth of current directory from main PST directory
pstVer =    'pstSETO';%   'pstV3p1'; % 'pstV2P3'; % 
pathParts = strsplit(pwd, filesep);
PSTpath = pathParts(1);

for pNdx = 2:max(size(pathParts))-folderDepth
    PSTpath = [char(PSTpath), filesep, char(pathParts(pNdx))];
end
PSTpath = [char(PSTpath), filesep, pstVer, filesep];

addpath(PSTpath)
save PSTpath.mat PSTpath pstVer caseName
clear folderDepth pathParts pNdx PSTpath

%% Run nonlinear simulation and store results
clear all; close all; clc
load PSTpath.mat

delete([PSTpath 'DataFile.m']); % ensure batch datafile is cleared
copyfile('datanp48.m',[PSTpath 'DataFile.m']); % copy system data file to batch run location

livePlotFlag = 1;
pssGainFix = 1;
s_simu_Batch %Run PST <- this is the main file to look at for simulation workings

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
save([pstVer, caseName, '.mat']); %Save simulation outputs

%% temp file clean up
delete('PSTpath.mat')
delete('sim_fle.mat')
