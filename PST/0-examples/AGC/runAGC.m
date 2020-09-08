% Test AGC two area case (AGC in 4 version only)

clear all; close all; clc
%% Add pst path to MATLAB
% generate relative path generically
folderDepth = 2; % depth of current directory from main PST directory
pstVer =  'PSTv4'; %'pstSETO'; % 'pstV2p3';%   'pstV3P1';%

% automatically handle global g usage
if strcmp(pstVer, 'pstSETO') || strcmp(pstVer, 'PSTv4')
    useGlobalG = true;
else
    useGlobalG = false;
end

pathParts = strsplit(pwd, filesep);
PSTpath = pathParts(1);

for pNdx = 2:max(size(pathParts))-folderDepth
    PSTpath = [char(PSTpath), filesep, char(pathParts(pNdx))];
end
PSTpath = [char(PSTpath), filesep, pstVer, filesep];

addpath(PSTpath)

save PSTpath.mat PSTpath pstVer useGlobalG
clear folderDepth pathParts pNdx PSTpath


%% Run nonlinear simulation and store results
clear all; close all; clc
load PSTpath.mat

% Test of variable time step
% available options:
% 'yes' - use Huen's method
% 'no' - use Variable time step methods described in d2a_AGC_VTS
% 'mixed' - use Variable time step methods described in d2a_AGC_mixed 
FTS = 'yes';

delete([PSTpath 'DataFile.m']); % ensure batch datafile is cleared

if strcmp(FTS,'no')
    copyfile('d2a_AGC_VTS.m',[PSTpath 'DataFile.m']); % copy system data file to batch run location
elseif strcmp(FTS,'mixed')
    copyfile('d2a_AGC_mixed.m',[PSTpath 'DataFile.m']); % copy system data file to batch run location
else
    copyfile('d2a_AGC.m',[PSTpath 'DataFile.m']); % copy system data file to batch run location
end

copyfile([PSTpath 'livePlot_2.m'],[PSTpath 'livePlot.m']); % specify plot operation
livePlotFlag = 1;

% Handle load modulation
if useGlobalG
    copyfile( 'ml_sig_smallStepG.m',[PSTpath 'ml_sig.m']); % For global G pstSETO
else
    error('AGC not included in PST version')
end

copyfile([PSTpath 'mac_sub_NEW2.m'],[PSTpath 'mac_sub.m']); % use updated model
copyfile([PSTpath 'pss2.m'],[PSTpath 'pss.m']); % use version 2 pss

s_simu %Run PST 

% reset original files
copyfile([PSTpath 'ml_sig_ORIG.m'],[PSTpath 'ml_sig.m']); % reset modulation file
copyfile([PSTpath 'livePlot_ORIG.m'],[PSTpath 'livePlot.m']); % reset live plot
copyfile([PSTpath 'mac_sub_ORIG.m'],[PSTpath 'mac_sub.m']); % use updated model
copyfile([PSTpath 'pss3.m'],[PSTpath 'pss.m']); % use version 3 pss

% delete PSTpath
delete PSTpath.mat

%% Save cleaned output data
caseName = 'AGC';
printFigs = 0;

if strcmp(FTS,'no')
    save([caseName,'testVTS2.mat']); %Save simulation outputs
elseif strcmp(FTS,'mixed')
    save([caseName,'testVTSmix.mat']); %Save simulation outputs
else
    save([caseName,'testFTS.mat']); %Save simulation outputs
end

if printFigs
    set(gcf,'color','w'); % to remove border of figure
    export_fig([caseName,'f1'],'-png'); % to print fig
end

%% Plotting ====================================================================
%% system frequencies
figure
plot(g.sys.t, g.area.area(1).aveF,'--')
hold on
plot(g.sys.t, g.area.area(2).aveF,'--')
plot(g.sys.t, g.sys.aveF,'k','linewidth',1.5)

grid on
title('Weighted Average Frequencies')
ylabel('Frequency [PU]')
xlabel('Time [sec]')
legend({'Area 1', 'Area 2', 'Total System'},'location','best')

if printFigs
    set(gcf,'color','w'); % to remove border of figure
    export_fig([caseName,'f2'],'-pdf'); % to print fig
end
%% AGC values
figure
% IC error
plot(g.sys.t, real(g.area.area(1).icA) - real(g.area.area(1).icA(1)),'.-','linewidth',2)
hold on
plot(g.sys.t, real(g.area.area(2).icA)- real(g.area.area(2).icA(1)),'.-','linewidth',2)

% RACE
plot(g.sys.t, g.agc.agc(1).race)
plot(g.sys.t, g.agc.agc(2).race)

% SACE
plot(g.sys.t, g.agc.agc(1).sace, ':','linewidth',2)
plot(g.sys.t, g.agc.agc(2).sace, ':','linewidth',2)

% ace2dist - signal actively being applied to area generators
plot(g.sys.t, g.agc.agc(1).ace2dist, '--' )
plot(g.sys.t, g.agc.agc(2).ace2dist, '--')

grid on
ylabel('MW [PU]')
xlabel('Time [sec]')
title('AGC values')
legend({'Area 1 IC error', 'Area 2 IC error', ...
    'Area 1 RACE', 'Area 2 RACE', ...
    'Area 1 SACE', 'Area 2 SACE', ...
    'Area 1 DACE', 'Area 2 DACE'} ...
    ,'location','best')

if printFigs
    set(gcf,'color','w'); % to remove border of figure
    export_fig([caseName,'f3'],'-pdf'); % to print fig
end

%% tg values
figure

plot(g.sys.t, g.tg.tg_sig(1,:))
hold on
plot(g.sys.t, g.tg.tg_sig(2,:))
plot(g.sys.t, g.tg.tg_sig(3,:))
plot(g.sys.t, g.tg.tg_sig(4,:))

grid on
ylabel('MW [PU]')
xlabel('Time [sec]')
title('Signal sent to governor Pref (tg\_sig)')
legend({'Gov 1', 'Gov 2','Gov 3', 'Gov 4'},'location','best')

if printFigs
    set(gcf,'color','w'); % to remove border of figure
    export_fig([caseName,'f4'],'-pdf'); % to print fig
end
%% change area generation
figure
plot(g.sys.t, real(g.area.area(1).totGen)-real(g.area.area(1).totGen(1)))
hold on
plot(g.sys.t, real(g.area.area(2).totGen)-real(g.area.area(2).totGen(1)))

grid on
ylabel('MW [PU]')
xlabel('Time [sec]')
title('Change in Area Generation')
legend({'Area 1', 'Area 2'},'location','best')

if printFigs
    set(gcf,'color','w'); % to remove border of figure
    export_fig([caseName,'f5'],'-pdf'); % to print fig
end

%% power to load bus via lmon
figure
plot(g.sys.t, real(g.lmon.line(1).sFrom)-real(g.lmon.line(1).sFrom(1)))
hold on
plot(g.sys.t, real(g.lmon.line(2).sFrom)-real(g.lmon.line(2).sFrom(1)))
plot(g.sys.t, imag(g.lmon.line(1).sFrom)-imag(g.lmon.line(1).sFrom(1)), '--')
plot(g.sys.t, imag(g.lmon.line(2).sFrom)-imag(g.lmon.line(2).sFrom(1)), '--')

grid on
ylabel('MW or MVAR [PU]')
xlabel('Time [sec]')
title('Change in Power Absorbed by Load Buses')
legend({'REAL - Bus 4', 'REAL - Bus 14','REAC - Bus 4', 'REAC - Bus 14'},'location','best')

if printFigs
    set(gcf,'color','w'); % to remove border of figure
    export_fig([caseName,'f6'],'-pdf'); % to print fig
end
%% monitored load bus voltages
figure
plot(g.sys.t, abs(g.bus.bus_v(3,:))-abs(g.bus.bus_v(3,1)))
hold on
plot(g.sys.t, abs(g.bus.bus_v(10,:))-abs(g.bus.bus_v(10,1)))

grid on
ylabel('Voltage [PU]')
xlabel('Time [sec]')
title('Change Load Bus Voltage')
legend({'Bus 4', 'Bus 14'},'location','best')

if printFigs
    set(gcf,'color','w'); % to remove border of figure
    export_fig([caseName,'f7'],'-pdf'); % to print fig
end


