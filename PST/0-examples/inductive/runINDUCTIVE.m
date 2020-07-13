% inductive machine/generator and inductive load case
% Tested as working in all Versions
% output the same in all versions
% commented out g.xxx variables required for SETO runs

clear all; close all; clc

%% Add pst path to MATLAB
% generate relative path generically
folderDepth = 2; % depth of current directory from main PST directory
pstVer =   'pstSETO'; % 'pstV2p3';%  'pstV3P1';%   
pathParts = strsplit(pwd, filesep);
PSTpath = pathParts(1);

for pNdx = 2:max(size(pathParts))-folderDepth
    PSTpath = [char(PSTpath), filesep, char(pathParts(pNdx))];
end
PSTpath = [char(PSTpath), filesep, pstVer, filesep];

addpath(PSTpath)
save PSTpath.mat PSTpath pstVer
clear folderDepth pathParts pNdx PSTpath

%% Run nonlinear simulation and store results
clear all; close all; clc
load PSTpath.mat
delete([PSTpath 'DataFile.m']); % ensure batch datafile is cleared
copyfile('data3mIg.m',[PSTpath 'DataFile.m']); % copy system data file to batch run location

% move inductive load
copyfile([PSTpath 'mac_ind2.m'],[PSTpath 'mac_ind.m']); % copy system data file to batch run location
% move modulation file
copyfile( 'ml_sig_smallStepG.m',[PSTpath 'ml_sig.m']); % For global G pstSETO
% copyfile( 'ml_sig_smallStep.m',[PSTpath 'ml_sig.m']); % for v 2.3 and 3.1


s_simu_Batch %Run PST <- this is the main file to look at for simulation workings

%reset inductive load file
copyfile([PSTpath 'mac_ind2.m'],[PSTpath 'mac_ind.m']); % copy system data file to batch run location
% reset modulation file
copyfile([PSTpath 'ml_sig_ORIG.m'],[PSTpath 'ml_sig.m']); % copy system data file to batch run location


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
save([pstVer,'INDnonLIN.mat']); %Save simulation outputs

%% PST linear system creation
clear all; close all;

svm_mgen_Batch

% MATLAB linear system creation using linearized PST results
tL = (0:0.001:10); % time to match PST d file time
modSig=zeros(1,size(tL,2)); % create blank mod signal same length as tL vector
modSig(find(tL>= 1 ))= 0.1; % mirror logic from exciterModSig into input vector
modSig(find(tL>= 2))= 0; % mirror logic from exciterModSig into input vector

bsys = b_lmod;
csys = [c_v;c_ang]; % inductive models output c is probably angle...

G = ss(a_mat,bsys,csys,zeros(size(csys,1),size(bsys,2))); % create system using pst matricies

y = lsim(G,modSig,tL); % run input into state space system

% collect bus voltage magnitudes and adjust by initial conditions
linV = y(:,1:size(c_v,1))'; % rotate into col vectors
linAng = y(:,size(c_v,1)+1:end)'; % collect and rotate angle data

% adjust data changes by initial conditions
for busN = 1:size(linV,1)
    linV(busN,:) = linV(busN,:) + bus_sol(busN,2);
    linAng(busN,:) = linAng(busN,:) + deg2rad(bus_sol(busN,3));
end

save linResults.mat tL linV linAng modSig 

%% plot comparisons
load PSTpath.mat
name = [pstVer,'INDnonLIN.mat'];
feval('load', name)
load linResults.mat

%% compare mod inputs
figure
hold on
plot(tL,modSig)
plot(t,g.lmod.lmod_sig,'--')
% plot(t,lmod_sig,'--')
legend('Linear','Non-Linear','location','best')
title('Governor Pref Modulation Signal')

%% compare bus voltage magnitude
figure
hold on
legNames={};
for busN=1:size(linV,1)
    plot(tL,linV(busN,:))
    legNames{end+1}= ['Bus ', int2str(busN), ' Linear'];
    plot(t,abs(g.sys.bus_v(busN,:)),'--')
%     plot(t,abs(bus_v(busN,:)),'--')
    legNames{end+1}= ['Bus ', int2str(busN), ' non-Linear'];
    
end
legend(legNames,'location','best')
title('Bus Voltage Magnitude')
xlabel('Time [sec]')
ylabel('Voltage [PU]')

%% compare bus angle
figure
hold on
legNames={};
for busN=1:size(linAng,1)
    plot(tL,linAng(busN,:))
    legNames{end+1}= ['Bus ', int2str(busN), ' Linear'];
    plot(t,angle(g.sys.bus_v(busN,:)),'--')
%     plot(t,angle(bus_v(busN,:)),'--')
    legNames{end+1}= ['Bus ', int2str(busN), ' non-Linear'];
    
end
legend(legNames,'location','best')
title('Bus Voltage Angle')
xlabel('Time [sec]')
ylabel('Angle [PU]')
