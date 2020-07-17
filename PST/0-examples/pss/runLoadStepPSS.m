% pss test
clear all; close all; clc

%% Add pst path to MATLAB
% generate relative path generically
folderDepth = 2; % depth of current directory from main PST directory
pstVer = 'pstSETO'; % 'pstV2p3';% 'pstV3P1';%   
pathParts = strsplit(pwd, filesep);
PSTpath = pathParts(1);

for pNdx = 2:max(size(pathParts))-folderDepth
    PSTpath = [char(PSTpath), filesep, char(pathParts(pNdx))];
end
PSTpath = [char(PSTpath), filesep, pstVer, filesep];

addpath(PSTpath)
save PSTpath.mat PSTpath
clear folderDepth pathParts pNdx PSTpath

% %% Select PSS model
% load PSTpath.mat
% delete([PSTpath 'pss.m']); % ensure batch datafile is cleared
% copyfile([PSTpath 'pss3.m'],[PSTpath 'pss.m']); 

%% Run nonlinear simulation and store results
clear all; close all; clc
load PSTpath.mat
delete([PSTpath 'DataFile.m']); % ensure batch datafile is cleared
copyfile('d_smallLoadStep.m',[PSTpath 'DataFile.m']); % copy system data file to batch run location

% Handle load modulation file placement etc...
delete([PSTpath 'ml_sig.m']); % ensure ml_sig file is empty
copyfile('ml_sig_smallStepG.m',[PSTpath 'ml_sig.m']); % copy simulation specific data file to batch run location

copyfile([PSTpath 'pss2.m'],[PSTpath 'pss.m']); % copy simulation specific data file to batch run location

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
        disp(ME.message)
        disp(vName)
    end

end
clear varNames vName zeroTest

%% Save cleaned output data
save('loadStepNONLIN.mat'); %Save simulation outputs

%% PST linear system creation
clear all; close all;

svm_mgen_Batch

%%  

%% MATLAB linear system creation using linearized PST results
tL = (0:0.001:20); % time to match PST d file time
modSig=zeros(1,size(tL,2)); % create blank mod signal same length as tL vector
modSig(find(tL> 1.0 ))= 0.01; % mirror logic from exciterModSig into input vector
modSig(find(tL>10))= 0; % mirror logic from exciterModSig into input vector

bsys = b_lmod;
csys = [c_v;c_spd;c_pm];
G = ss(a_mat,bsys,csys,zeros(size(csys,1),size(bsys,2))); % create system using pst matricies

y = lsim(G,modSig,tL); % run input into state space system

% collect bus voltage magnitudes and adjust by initial conditions
linV = y(:,1:4)'; % rotate into col vectors
for busN = 1:size(linV,1)
    linV(busN,:) = linV(busN,:) + bus_sol(busN,2);
end

% collect machine speeds and adjust by initial condition
linSpd = y(:,5:6)'+ 1.0; % rotate into col vectors

% collect pm...
linPm = y(:,7:8)';% rotate to vector
linPm(1,:)= linPm(1,:)+ g.mac.pmech(1,1);
linPm(2,:)= linPm(2,:)+ g.mac.pmech(2,1);


save linResults.mat tL linV linSpd modSig linPm

%% Clean up load modulation file alterations...
load PSTpath
delete([PSTpath 'ml_sig.m']); % remove simulation specific ml_sig file
copyfile([PSTpath 'ml_sig_ORIG.m'],[PSTpath 'ml_sig.m']); % Replace original file
clear all

%% temp file clean up
delete('PSTpath.mat')
delete('sim_fle.mat')

%% plot comparisons
load loadStepNONLIN.mat
load linResults.mat

%% compare mod inputs
figure
hold on
plot(tL,modSig)
plot(g.sys.t,g.lmod.lmod_sig,'--')
%plot(t,pm_sig,'--')
legend('Linear','Non-Linear','location','best')
title('Governor Pref Modulation Signal')


%% compare machine speeds
figure
hold on
legNames={};
for busN=1:size(linSpd,1)
    plot(tL,linSpd(busN,:))
    legNames{end+1}= ['Gen Speed ', int2str(busN), ' Linear'];
    plot(g.sys.t,g.mac.mac_spd(busN,:),'--')
    legNames{end+1}= ['Gen Speed ', int2str(busN), ' non-Linear'];
    
end
legend(legNames,'location','best')
title('Machine Speeds')
xlabel('Time [sec]')
ylabel('Speed [PU]')

%% compare machine power
figure
hold on
legNames={};
for busN=1:size(linPm,1)
    plot(tL,linPm(busN,:))
    legNames{end+1}= ['Gen Pm ', int2str(busN), ' Linear'];
    plot(g.sys.t,g.mac.pmech(busN,:),'--')
    legNames{end+1}= ['Gen Pm ', int2str(busN), ' non-Linear'];
end
legend(legNames,'location','best')
title('Machine Mechanical Power')
xlabel('Time [sec]')
ylabel('Mechanical Power [PU MW]')

%% compare bus voltage magnitude
figure
hold on
legNames={};
for busN=1:size(linV,1)
    plot(tL,linV(busN,:))
    legNames{end+1}= ['Bus ', int2str(busN), ' Linear'];
    plot(g.sys.t,abs(g.sys.bus_v(busN,:)),'--')
    legNames{end+1}= ['Bus ', int2str(busN), ' non-Linear'];
    
end
legend(legNames,'location','best')
title('Bus Voltage Magnitude')
xlabel('Time [sec]')
ylabel('Voltage [PU]')
