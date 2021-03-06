% Example of miniwecc test
% pstV2P3 - seems to work - 50.5244 (59.3868s) L
% pstSETO - seems to work - 18.9074 (37.0102s) L
% pstV3p1 - seems to work - 50.2106 (58.7085s) L
% PSTv3 - has no batch run.... assumed the same as v3p1.
%

clear all; close all; clc
caseName = 'PSS';
scenario = 'C';% L line, F colstrip fault, C cascade?

%% Add pst path to MATLAB
% generate relative path generically
folderDepth = 3; % depth of current directory from main PST directory
pstVer =    'PSTv4';  %'pstSETO';   % 'pstV3p1'; %  'pstV2P3'; %
pathParts = strsplit(pwd, filesep);
PSTpath = pathParts(1);

for pNdx = 2:max(size(pathParts))-folderDepth
    PSTpath = [char(PSTpath), filesep, char(pathParts(pNdx))];
end
PSTpath = [char(PSTpath), filesep, pstVer, filesep];

addpath(PSTpath)
save PSTpath.mat PSTpath pstVer caseName scenario
clear folderDepth pathParts pNdx PSTpath scenario

%% Run nonlinear simulation and store results
clear all; close all; clc
load PSTpath.mat
delete([PSTpath 'DataFile.m']); % ensure batch datafile is cleared

switch scenario
    case 'L'
        copyfile('d_miniWECC_V3C_C3_6_C_AlbertaSw_LineOpen_CurrentLoads.m',[PSTpath 'DataFile.m']); % copy system data file to batch run location
    case 'F'
        copyfile('d_minniWECC_ColstripFault.m',[PSTpath 'DataFile.m']); % copy system data file to batch run location
    case 'C'
        copyfile('d_minniWECC_CascadedFault5.m',[PSTpath 'DataFile.m']); % copy system data file to batch run location
    otherwise
        error('invalid scenario')
end

%
% % Handle load modulation file placement etc...
% copyfile([PSTpath 'ml_sig.m'],[PSTpath 'ml_sig_ORIG_TMP.m']); % save copy of original ml_sig file
% delete([PSTpath 'ml_sig.m']); % ensure ml_sig file is empty
% copyfile('ml_sig_loadStep.m',[PSTpath 'ml_sig.m']); % copy simulation specific data file to batch run location
%
% % Handle turbine modulation file placement etc...
% copyfile([PSTpath 'mtg_sig.m'],[PSTpath 'mtg_sig_ORIG_TMP.m']); % save copy of original ml_sig file
% delete([PSTpath 'mtg_sig.m']); % ensure ml_sig file is empty
% copyfile('mtg_sig_PrefStep.m',[PSTpath 'mtg_sig.m']); % copy simulation specific data file to batch run location

livePlotFlag = 1;
pssGainFix = 0;

if strcmp(pstVer , 'PSTv4')
    s_simu
else
    s_simu_Batch %Run PST 
end
%% Clean up modulation file alterations.
% delete([PSTpath 'ml_sig.m']); % remove simulation specific ml_sig file
% copyfile([PSTpath 'ml_sig_ORIG_TMP.m'],[PSTpath 'ml_sig.m']); % Replace original file
% delete([PSTpath 'ml_sig_ORIG_TMP.m']); % delete temporary file
%
% % turbin mod
% delete([PSTpath 'mtg_sig.m']); % remove simulation specific ml_sig file
% copyfile([PSTpath 'mtg_sig_ORIG_TMP.m'],[PSTpath 'mtg_sig.m']); % Replace original file
% delete([PSTpath 'mtg_sig_ORIG_TMP.m']); % delete temporary file



%% Save cleaned output data
save([pstVer, caseName, scenario, '.mat']); %Save simulation outputs

%% temp file clean up
delete('PSTpath.mat')
