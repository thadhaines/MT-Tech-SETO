% Example of one machine infinite bus line trip using d_OneMacInfBusXX
% smaller system used to more easily study inner workings of PST 
%% 04
% Experimentation with load modulation
% Load step on bus 2 of +0.25 PU at t=1 where lmod TS = 0.1 second
% One steam gov gen, One hydro gov.
% Hydro gov twice size as Steam.

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
save PSTpath.mat PSTpath
clear folderDepth pathParts pNdx PSTpath

%% Run nonlinear simulation and store results
clear all; close all; clc
load PSTpath.mat
delete([PSTpath 'DataFile.m']); % ensure batch datafile is cleared
copyfile('d_OneMacInfBus04.m',[PSTpath 'DataFile.m']); % copy system data file to batch run location

% Handle load modulation file placement etc...
copyfile([PSTpath 'ml_sig.m'],[PSTpath 'ml_sig_ORIG_TMP.m']); % save copy of original ml_sig file
delete([PSTpath 'ml_sig.m']); % ensure ml_sig file is empty
copyfile('ml_sig_OMIB3.m',[PSTpath 'ml_sig.m']); % copy simulation specific data file to batch run location

s_simu_Batch %Run PST <- this is the main file to look at for simulation workings

%% Clean up load modulation file alterations...
delete([PSTpath 'ml_sig.m']); % remove simulation specific ml_sig file
copyfile([PSTpath 'ml_sig_ORIG_TMP.m'],[PSTpath 'ml_sig.m']); % Replace original file
delete([PSTpath 'ml_sig_ORIG_TMP.m']); % delete temporary file

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
        disp(ME.message)
        disp(vName)
    end

end
clear varNames vName zeroTest

%% Save cleaned output data
save('OneMacInfBus04.mat'); %Save simulation outputs

%% temp file clean up
delete('PSTpath.mat')
delete('sim_fle.mat')

%% Plotting init and function call
clear all % ensure only saved data is plotted
load('OneMacInfBus04.mat')

%pstMegaPlot
disp('fin')