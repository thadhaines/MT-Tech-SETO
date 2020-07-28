%% File to New England 39 bus system
clear all; close all; clc
caseName = 'datane_hiskens';

%% Add pst path to MATLAB
% generate relative path generically
folderDepth = 3; % depth of current directory from main PST directory
pstVer =   'pstSETO';%  'pstV3p1'; % 'pstV2P3'; %  
pathParts = strsplit(pwd, filesep);
PSTpath = pathParts(1);

for pNdx = 2:max(size(pathParts))-folderDepth
    PSTpath = [char(PSTpath), filesep, char(pathParts(pNdx))];
end
PSTpath = [char(PSTpath), filesep, pstVer, filesep];

addpath(PSTpath)
addpath([PSTpath, 'test', filesep]) % to handle new functionized code
save PSTpath.mat PSTpath pstVer caseName
clear folderDepth pathParts pNdx PSTpath

%% Run nonlinear simulation and store results
clear all; close all; clc
load PSTpath.mat

delete([PSTpath 'DataFile.m']); % ensure batch datafile is cleared
copyfile([caseName, '.m'],[PSTpath 'DataFile.m']); % copy system data file to batch run location

copyfile([PSTpath 'pss3.m'],[PSTpath 'pss.m']); % specify pss model
copyfile([PSTpath 'mac_sub_NEW2.m'],[PSTpath 'mac_sub.m']); % specify machine model
copyfile([PSTpath 'livePlot_1.m'],[PSTpath 'livePlot.m']); % specify plot operation
livePlotFlag = 1;

% 15s - slower during transients - faster when no action.
% 113 - works well during transierts, slower during no action 
% ode23s - many iterations per step (~200) not viable
% ode23tb - occasionally hundereds of iterations, sometimes not... decent
% ode23 - similar to 23tb, timstep doesn't get very large
% ode23t - works...

odeName = 'ode23t';

% s_simu_Batch %Run PST with original format
% s_simu_BatchTestF %Run PST functionalized test
s_simu_BatchVTS %Run PST with variable timestep

copyfile([PSTpath 'pss3.m'],[PSTpath 'pss.m']); % reset pss
copyfile([PSTpath 'mac_sub_ORIG.m'],[PSTpath 'mac_sub.m']); % reset mac_sub
copyfile([PSTpath 'livePlot_ORIG.m'],[PSTpath 'livePlot.m']); % reset live plot


%% Save cleaned output data
save(['VTS',pstVer, caseName, '.mat']); %Save simulation outputs
save(['VTS',odeName, '.mat']); %Save simulation outputs

%% temp file clean up
delete('PSTpath.mat')

compareVTSandFTS
