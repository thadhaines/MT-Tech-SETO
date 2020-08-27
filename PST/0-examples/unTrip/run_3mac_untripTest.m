% Run 3mac_untripTest
clear all; close all; clc

%% Add correct PST verstion path to MATLAB in a generic way
folderDepth = 2; % depth of current directory from main PST directory
pstVer =   'PSTv4'; % 'pstV3p1'; % 'pstV2P3'; % 'pstSETO';%
pathParts = strsplit(pwd, filesep);
PSTpath = pathParts(1);

for pNdx = 2:max(size(pathParts))-folderDepth
    PSTpath = [char(PSTpath), filesep, char(pathParts(pNdx))];
end
PSTpath = [char(PSTpath), filesep, pstVer, filesep];

addpath(PSTpath)
save PSTpath.mat PSTpath pstVer
clear folderDepth pathParts pNdx PSTpath scenario

%% Run nonlinear simulation and store results
clear all; close all; clc
load PSTpath.mat

copyfile('d_3mac_untripTest.m',[PSTpath 'DataFile.m']); % system data file
copyfile([PSTpath 'livePlot_1.m'],[PSTpath 'livePlot.m']); % set live plot

% % Move correct modulation
% if strcmp(pstVer, 'pstSETO') || strcmp(pstVer, 'PSTv4')
%     copyfile('ml_sig_smallStepG.m', [PSTpath 'ml_sig.m']);
% else
%     error('not handled yet')
% end

% Move correct trip logic files
if strcmp(pstVer, 'pstSETO') || strcmp(pstVer, 'PSTv4')
    copyfile('mac_trip_logic_Gen_3_G.m', [PSTpath 'mac_trip_logic.m']);
else
    error('not handled yet')
end

livePlotFlag = 1; % not ideal for extended term sim - may cause crashes/freezing

if strcmp(pstVer , 'PSTv4')
    s_simu
else
    s_simu_Batch %Run PST 
end

%% Clean up modulation file alterations.
copyfile([PSTpath 'mac_trip_logic_ORIG.m'], [PSTpath 'mac_trip_logic.m']);
copyfile([PSTpath 'livePlot_ORIG.m'],[PSTpath 'livePlot.m']);
copyfile([PSTpath 'ml_sig_ORIG.m'],[PSTpath 'ml_sig.m']);

%% Save cleaned output data
save('3mac_untripTest.mat'); %Save simulation outputs

%% temp file clean up
delete('PSTpath.mat')

%% =============================================================== 

%% Plots

plotCell = { ...
    %f1, f2
    'mac','mac_spd';
    'mac','pelect';    
    'mac','qelect';  
    'mac','cur_re'; 
    'mac','cur_im';
    'mac','mac_ang';
    };

% nS = find(g.sys.t > 24);
% nS = nS(1);
% nE = find(g.sys.t > 28);
% nE = nE(1);
axLim = [24,28];

for n = 1:size(plotCell,1)
    f1 = plotCell{n,1};
    f2 = plotCell{n,2};
figure
plot(g.sys.t, g.(f1).(f2))
xlabel('Time [sec]')
%xlim(axLim)
title([f1,'.',f2], 'Interpreter','None')
end

% mac trip flag held... probably has to do with a latching thing...

