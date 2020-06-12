% Example of two machine system with exciter modulation
% smaller system used to more easily study inner workings of PST
% One steam gov gen, One hydro gov.
% only 1 exciter
% Hydro gov twice size as Steam.

% NOTE: live plotting disabled in d_ file to test speed up and better compare MATLAB to Octave performance
% (Octave seems a bit faster when plotting is removed from simulation run)
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
copyfile('d_OMIB_EXC.m',[PSTpath 'DataFile.m']); % copy system data file to batch run location

% Handle exciter modulation file placement etc.
copyfile([PSTpath 'mexc_sig.m'],[PSTpath 'mexc_sig_ORIG_TMP.m']); % save copy of original ml_sig file
delete([PSTpath 'mexc_sig.m']); % ensure ml_sig file is empty
copyfile('exciterModSig.m',[PSTpath 'mexc_sig.m']); % copy simulation specific data file to batch run location

s_simu_Batch %Run PST non-linear sim

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
save('exciterTest.mat'); %Save simulation outputs

%% PST linear system creation
clear all; close all;
svm_mgen_Batch

%% MATLAB linear system creation using linearized PST results
tL = (0:0.01:20); % time to match PST d file time
excSig=zeros(1,size(tL,2)); % create blank mod signal same length as tL vector
excSig(find(tL>0.5))= 0.01; % mirror logic from exciterModSig into input vector
G = ss(a_mat,b_vr,[c_v;c_spd],zeros(6,1)); % create system using pst matricies

y = lsim(G,excSig,tL); % run input into state space system

% collect bus voltage magnitudes and adjust by initial conditions
linV = y(:,1:4)'; % rotate into col vectors
for busN = 1:size(linV,1)
    linV(busN,:) = linV(busN,:) + bus_sol(busN,2);
end

% collect machine speeds and adjust by initial condition
linSpd = y(:,5:6)'+1.0; % rotate into col vectors
save linResults.mat tL linV linSpd excSig
clear all

%% Clean up modulation file alterations...
load PSTpath.mat
delete([PSTpath 'mexc_sig.m']); % remove simulation specific ml_sig file
copyfile([PSTpath 'mexc_sig_ORIG_TMP.m'],[PSTpath 'mexc_sig.m']); % Replace original file
delete([PSTpath 'mexc_sig_ORIG_TMP.m']); % delete temporary file

%% temp file clean up
delete('PSTpath.mat')
delete('sim_fle.mat')

%% plot comparisons
load exciterTest.mat
load linResults.mat
%% compare Inputs to exciter
figure
hold on
plot(tL,excSig)
plot(t,exc_sig,'--')
legend('Linear','Non-Linear','location','best')
title('Exciter Modulation Signal')
%% compare bus voltage magnitude
figure
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
%% compare machine speeds
figure
hold on
legNames={};
for busN=1:size(linSpd,1)
    plot(tL,linSpd(busN,:))
    legNames{end+1}= ['Gen Speed ', int2str(busN), ' Linear'];
    plot(t,mac_spd(busN,:),'--')
    legNames{end+1}= ['Gen Speed ', int2str(busN), ' non-Linear'];
    
end
legend(legNames,'location','best')
title('Machine Speeds')