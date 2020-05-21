% Example of one machine infinite bus line trip using d_OneMacInfBus
% smaller system used to more easily study inner workings of PST 
% Include work arounds for Octave by running octaveComp (saved in PST main directory)
% tested as working with octave 5.2.0 on 5/21/20

clear all; close all; clc

%% Add pst path to current directory
% generate relative path generically
folderDepth = 4; % depth of current directory from main PST directory

pathParts = strsplit(pwd, filesep);
PSTpath = pathParts(1);

for pNdx = 2:max(size(pathParts))-folderDepth
    PSTpath = [char(PSTpath), filesep, char(pathParts(pNdx))];
end
PSTpath = [char(PSTpath), filesep];

addpath(PSTpath)
save PSTpath PSTpath
clear folderDepth pathParts pNdx PSTpath

%% Run nonlinear simulation and store results
clear all; close all; clc
load PSTpath
delete([PSTpath 'DataFile.m']); % ensure batch datafile is cleared
copyfile('d_OneMacInfBus.m',[PSTpath 'DataFile.m']); % copy system data file to batch run location
s_simu_Batch %Run PST <- this is the main file to look at for simulation workings

%% Simulation variable cleanup
% Clear any varables that contain only zeros
varNames = who()'; % all variable names in workspace
clearedVars = {}; % cell to hold names of deleted 'all zero' variables

for vName = varNames
    zeroTest = eval(sprintf('all(%s(:)==0)', vName{1})); % check if all zeros
    if zeroTest
        eval(sprintf('clear %s',vName{1}) ); % clear variable
        clearedVars{end+1} = vName{1}; % add name to cell for reference
    end
end
clear varNames vName zeroTest

%% Save cleaned output data
save('OneMacInfBus01Oct.mat'); %Save simulation outputs

%% temp file clean up
delete('PSTpath.mat')
delete('sim_fle.mat')

%% Plotting init and function call
clear all % ensure only saved data is plotted
load('OneMacInfBus01Oct.mat')

pstMegaPlot