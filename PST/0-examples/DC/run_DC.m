% DC line load case
% Tested as 'working' in all Versions
% Linear simulation doesn't match non-linear output - probably user error

clear all; close all; clc
%% Add pst path to MATLAB
% generate relative path generically
folderDepth = 2; % depth of current directory from main PST directory
pstVer =  'PSTv4'; %   'pstSETO'; % 'pstV2p3';%  'pstV3P1';% 

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
addpath([PSTpath, 'test', filesep]) % to handle new functionized code
save PSTpath.mat PSTpath pstVer useGlobalG
clear folderDepth pathParts pNdx PSTpath

%% Run nonlinear simulation and store results
clear all; close all; clc
load PSTpath.mat
delete([PSTpath 'DataFile.m']); % ensure batch datafile is cleared
copyfile('d_testdcREF.m',[PSTpath 'DataFile.m']); % copy system data file to batch run location

% Handle load modulation
if useGlobalG
    copyfile( 'ml_sig_smallStepG.m',[PSTpath 'ml_sig.m']); % For global G pstSETO
else
    copyfile( 'ml_sig_smallStep.m',[PSTpath 'ml_sig.m']); % for v 2.3 and 3.1
end

if strcmp(pstVer, 'pstSETO') || strcmp(pstVer, 'PSTv4')
    s_simu
else
    s_simu_Batch %Run PST
end
% % reset modulation file
copyfile([PSTpath 'ml_sig_ORIG.m'],[PSTpath 'ml_sig.m']);


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
save([pstVer,'DCnonLIN.mat']); %Save simulation outputs

%% PST linear system creation
clear all; close all;

svm_mgen_Batch

% MATLAB linear system creation using linearized PST results
tL = (0:0.001:5); % time to match PST d file time
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
load PSTpath
for busN = 1:size(linV,1)
    if strcmp(pstVer, 'pstSETO') || strcmp(pstVer, 'PSTv4')
        linV(busN,:) = linV(busN,:) + g.bus.bus(busN,2);
        linAng(busN,:) = linAng(busN,:) + deg2rad(g.bus.bus(busN,3));
    else
        linV(busN,:) = linV(busN,:) + bus(busN,2);
        linAng(busN,:) = linAng(busN,:) + deg2rad(bus(busN,3));
    end
    
end

load PSTpath.mat
save([pstVer,'linResults.mat'], 'tL', 'linV', 'linAng', 'modSig')

%% plot comparisons
name = [pstVer,'DCnonLIN.mat'];
feval('load', name)
load([pstVer,'linResults.mat'])

%% compare mod inputs
figure
hold on
plot(tL,modSig)

if useGlobalG
    plot(g.sys.t,g.lmod.lmod_sig,'--')
else
    plot(t,lmod_sig,'--')
end

legend('Linear','Non-Linear','location','best')
title('Governor Pref Modulation Signal')

%% compare bus voltage magnitude
figure
hold on
legNames={};
for busN=1:size(linV,1)
    plot(tL,linV(busN,:))
    legNames{end+1}= ['Bus ', int2str(busN), ' Linear'];
    
    if useGlobalG
        plot(g.sys.t,abs(g.bus.bus_v(busN,:)),'--')
    else
        plot(t,abs(bus_v(busN,:)),'--')
    end
    
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
    
    if useGlobalG
        plot(g.sys.t,angle(g.bus.bus_v(busN,:)),'--')
    else
        plot(t,angle(bus_v(busN,:)),'--')
    end
    
    legNames{end+1}= ['Bus ', int2str(busN), ' non-Linear'];
    
end
legend(legNames,'location','best')
title('Bus Voltage Angle')
xlabel('Time [sec]')
ylabel('Angle [PU]')
