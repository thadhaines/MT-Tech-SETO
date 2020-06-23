%% runBatchTests.m
%   Script to run PST tests aimed at testing specific models
%   Modified from full example cases with shortened simulation times.

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

%% create .mat for completed tests assuming workspace repeatidly cleared
compTest = {};
save compTest.mat compTest

%% Run exciter type 0 test
runEXCtype0

load compTest.mat
% append Test name and data prefix to running completed test Cell
compTest{length(compTest)+1} = {'Exciter - Type 0 - Simple Exciter - smpexc', 'excT0', 'Exciter Modulation [V PU]'};
save compTest.mat compTest

%% Run exciter type 1 test
runEXCtype1

load compTest.mat
% append Test name and data prefix to running completed test Cell
compTest{length(compTest)+1} = {'Exciter - Type 1 - DC Exciter Type 1 - exc_dc12', 'excT1', 'Exciter Modulation [V PU]'};
save compTest.mat compTest

%% Run exciter type 2 test
runEXCtype2

load compTest.mat
% append Test name and data prefix to running completed test Cell
compTest{length(compTest)+1} = {'Exciter - Type 2 - DC Exciter Type 2 - exc_dc12', 'excT2', 'Exciter Modulation [V PU]'};
save compTest.mat compTest

%% Run exciter type 3 test
runEXCtype3

load compTest.mat
% append Test name and data prefix to running completed test Cell
compTest{length(compTest)+1} = {'Exciter - Type 3 - ST3 model - exc_st3', 'excT3', 'Exciter Modulation [V PU]'};
save compTest.mat compTest

%% Run exciter type 4 test
runEXCtype4

load compTest.mat
% append Test name and data prefix to running completed test Cell
compTest{length(compTest)+1} = {'Exciter - Type 4 - Simple exciter with PI - smppi', 'excT4', 'Exciter Modulation [V PU]'};
save compTest.mat compTest


%% Create plots for each test
clear all; close all; clc
load compTest.mat
for ut=1:length(compTest)
    % load data
    caseName = compTest{ut}{1};
    dataPrefix = compTest{ut}{2};
    modSigName = compTest{ut}{3};
    linDataN = [ dataPrefix, 'LIN.mat'];
    nlDataN = [dataPrefix, 'NL.mat'];
    feval('load',linDataN)
    feval('load', nlDataN)
    
    %% compare modulation signal
    figure('Name',caseName)
    subplot(2, 2, 1)
    hold on
    plot(tL,modSig)
    plot(t,modSigNL,'--')
    legend('Linear','Non-Linear','location','best')
    title('Modulation Signal')
    xlabel('Time [sec]')
    ylabel(modSigName)
    
    %% compare bus voltage magnitude
    subplot(2, 2, 2)
    hold on
    legNames={};
    for busN=1:size(linV,1)
        plot(tL,linV(busN,:))
        legNames{end+1}= ['Bus ', int2str(busN), ' Linear'];
        plot(t,abs(bus_v(busN,:)),'--')
        legNames{end+1}= ['Bus ', int2str(busN), ' non-Linear'];
        
    end
    legend(legNames,'location','best')
    title('Bus Voltage Magnitude')
    xlabel('Time [sec]')
    ylabel('Voltage [PU]')
    
    %% compare machine speeds
    subplot(2, 2, 3)
    hold on
    legNames={};
    for busN=1:size(linSpd,1)
        plot(tL,linSpd(busN,:))
        legNames{end+1}= ['Gen Speed ', int2str(busN), ' Linear'];
        plot(t,g.mac.mac_spd(busN,:),'--')
        legNames{end+1}= ['Gen Speed ', int2str(busN), ' non-Linear'];
        
    end
    legend(legNames,'location','best')
    title('Machine Speeds')
    xlabel('Time [sec]')
    ylabel('Speed [PU]')
    
    %% compare machine power
    subplot(2, 2, 4)
    hold on
    legNames={};
    for busN=1:size(linPm,1)
        plot(tL,linPm(busN,:))
        legNames{end+1}= ['Gen Pm ', int2str(busN), ' Linear'];
        plot(t,g.mac.pmech(busN,:),'--')
        legNames{end+1}= ['Gen Pm ', int2str(busN), ' non-Linear'];
    end
    legend(legNames,'location','best')
    title('Machine Mechanical Power')
    xlabel('Time [sec]')
    ylabel('Mechanical Power [PU MW]')
    
    %% clear data
    feval('clear',linDataN, nlDataN)
end
%% Clean up Temp .mat files
delete compTest.mat PSTpath.mat sim_fle.mat