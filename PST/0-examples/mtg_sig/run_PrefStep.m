% Example of two machine load step with govs
% smaller system used to more easily study inner workings of PST

% Experimentation with gov input (Pref) modulation

clear all; close all; clc

%% Add pst path to MATLAB
% generate relative path generically
folderDepth = 2; % depth of current directory from main PST directory
pstVer = 'PSTv4'; % 'pstSETO'; % 'pstV3p1';% 'pstV2P3';   %
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
copyfile('d_PrefStep.m',[PSTpath 'DataFile.m']); % copy system data file to batch run location

% Handle turbine modulation file placement etc...
if strcmp('PSTv4', pstVer) || strcmp('pstSETO', pstVer)
    copyfile('mtg_sig_PrefStepG.m',[PSTpath 'mtg_sig.m']);
else
    copyfile('mtg_sig_PrefStep.m',[PSTpath 'mtg_sig.m']);
end

% run non-linear simulation
if strcmp('PSTv4', pstVer)
    s_simu
else
    s_simu_Batch
end

%% Save output data
save('prefStepNonLin.mat'); %Save simulation outputs

%% PST linear system creation
clear all; close all;
svm_mgen_Batch

%% MATLAB linear system creation using linearized PST results
tL = (0:0.01:15); % time to match PST d file time
modSig=zeros(1,size(tL,2));         % create blank mod signal same length as tL vector
modSig(find(tL>0.5))= 0.05;         % mirror logic from input vector
modSig(find(tL>5))= -0.05;          % mirror logic from input vector
modSig(find(tL>10))= 0;             % mirror logic from input vector
G = ss(a_mat,b_pr,[c_v;c_spd],zeros(6,1)); % create system using pst matricies

y = lsim(G,modSig,tL); % run input into state space system

% collect bus voltage magnitudes and adjust by initial conditions
linV = y(:,1:4)'; % rotate into col vectors

load PSTpath
for busN = 1:size(linV,1)
    if strcmp('PSTv4', pstVer) || strcmp('pstSETO', pstVer)
        linV(busN,:) = linV(busN,:) + g.bus.bus(busN,2);
    else
        linV(busN,:) = linV(busN,:) + bus(busN,2);
    end
end

% collect machine speeds and adjust by initial condition
linSpd = y(:,5:6)'+ 1.0; % rotate into col vectors
save linResults.mat tL linV linSpd modSig
clear all
%% Clean up modulation file alterations.
% turbine governor pref mod
load PSTpath.mat
copyfile([PSTpath 'mtg_sig_ORIG.m'],[PSTpath 'mtg_sig.m']); % Replace original file

%% load data for plot comparisons
load linResults.mat
load prefStepNonLin.mat

%% compare Inputs to exciter
figure
hold on
plot(tL,modSig)
if strcmp('PSTv4', pstVer) || strcmp('pstSETO', pstVer)
    plot(g.sys.t,g.tg.tg_sig,'--')
else
    plot(t, tg_sig, '--')
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
    if strcmp('PSTv4', pstVer) || strcmp('pstSETO', pstVer)
        plot(g.sys.t,abs(g.bus.bus_v(busN,:)),'--')
    else
        plot(t,abs(bus_v(busN,:)),'--')
    end
    
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
    if strcmp('PSTv4', pstVer) || strcmp('pstSETO', pstVer)
        plot(g.sys.t, g.mac.mac_spd(busN,:),'--')
    else
        plot(t, mac_spd(busN,:),'--')
    end
    
    legNames{end+1}= ['Gen Speed ', int2str(busN), ' non-Linear'];
    
end
legend(legNames,'location','best')
title('Machine Speeds')

%% temp file clean up
delete('PSTpath.mat')
