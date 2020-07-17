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
save PSTpath.mat PSTpath pstVer caseName
clear folderDepth pathParts pNdx PSTpath

%% Run nonlinear simulation and store results
clear all; close all; clc
load PSTpath.mat

delete([PSTpath 'DataFile.m']); % ensure batch datafile is cleared
copyfile([caseName, '.m'],[PSTpath 'DataFile.m']); % copy system data file to batch run location

copyfile([PSTpath 'pss3.m'],[PSTpath 'pss.m']); % copy system data file to batch run location
copyfile([PSTpath 'mac_sub_NEW2.m'],[PSTpath 'mac_sub.m']); % copy system data file to batch run location
livePlotFlag = 1;

s_simu_Batch %Run PST <- this is the main file to look at for simulation workings

copyfile([PSTpath 'pss2.m'],[PSTpath 'pss.m']); % copy system data file to batch run location
copyfile([PSTpath 'mac_sub_ORIG.m'],[PSTpath 'mac_sub.m']); % copy system data file to batch run location

% %% Simulation variable cleanup
% % Clear any varables that contain only zeros
% varNames = who()'; % all variable names in workspace
% clearedVars = {}; % cell to hold names of deleted 'all zero' variables
% 
% for vName = varNames
%     try
%         zeroTest = eval(sprintf('all(%s(:)==0)', vName{1})); % check if all zeros
%         if zeroTest
%             eval(sprintf('clear %s',vName{1}) ); % clear variable
%             clearedVars{end+1} = vName{1}; % add name to cell for reference
%         end
%     catch ME
%         % gets called for structs... (global g)
%         disp(ME.message)
%         disp(vName)
%     end
% end
% clear varNames vName zeroTest

%% Save cleaned output data
save([pstVer, caseName, '.mat']); %Save simulation outputs

%% temp file clean up
delete('PSTpath.mat')
delete('sim_fle.mat')

%% plots to match dynamic response in hiskens pdf
% assumes pstSETO
% Generator Angle - ref gen 1 as per pdf
figure
for gen=1:g.mac.n_mac
        switch gen
        case 1
            yLIM = [-1, 1];
        case {2, 3, 4, 6, 7, 8, 9}
            yLIM = [-1,3];
        case 5
            yLIM = [-2, 4];
        case 10
            yLIM = [-1, 2];
        end
    subplot(5,2,gen)
    plot(g.sys.t, g.mac.mac_ang(gen,:)-g.mac.mac_ang(1,:),'k', 'linewidth',1.25)
    titleSTR = ['Gen',int2str(gen)];
    title(titleSTR)
    grid on
    ylim(yLIM)
    ylabel('Delta [rad]')
    xlabel('Time [sec]')
end

% plot to match machine response
figure
for gen=1:g.mac.n_mac
        switch gen
        case {1 , 2 , 9 , 10}
            yLIM = [-0.5,1.5];
        case {3 , 8 }
            yLIM = [-0.5, 1];
        case { 4 , 6 , 7 }
            yLIM = [-1, 2];
        case 5
            yLIM = [-2, 4];
        end
    subplot(5,2,gen)
    plot(g.sys.t, (g.mac.mac_spd(gen,:)-1)*100,'k', 'linewidth',1.25 )
    titleSTR = ['Gen',int2str(gen)];
    title(titleSTR)
    ylim(yLIM)
    grid on
    ylabel('Omega [%PU]')
    xlabel('Time [sec]')
end
clear yLIM titleSTR gen % variables with script plot

clear PSTpath caseName % script variables

clear pssV3 gov_flag % set in d file
