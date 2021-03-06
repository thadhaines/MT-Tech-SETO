% Example of two machine Pmech step with govs
% smaller system used to more easily study inner workings of PST

clear all; close all; clc

%% Add pst path to MATLAB
% generate relative path generically
folderDepth = 2; % depth of current directory from main PST directory
pstVer =  'PSTv4'; % 'pstSETO'; % 'pstV2P3'; %  'pstV3p1'; %
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
delete([PSTpath 'DataFile.m']);                     % ensure batch datafile is cleared
copyfile('d_2machineGov.m',[PSTpath 'DataFile.m']); % copy system data file to batch run location

% Handle modulation file placement etc...
if strcmp('PSTv4', pstVer) || strcmp('pstSETO', pstVer)
    copyfile('mpm_sig_PmStepG.m',[PSTpath 'mpm_sig.m']);
else
    copyfile('mpm_sig_PmStep.m',[PSTpath 'mpm_sig.m']);
end

% Use 'new' mac_sub model
copyfile([PSTpath 'mac_sub_NEW.m'],[PSTpath 'mac_sub.m']); % copy simulation specific data file to batch run location

if strcmp('PSTv4', pstVer)
    s_simu
else
    s_simu_Batch
end

%% Save cleaned output data
save('pMechStepNonLin.mat'); %Save simulation outputs

%% PST linear system creation
clear all; close all;
svm_mgen_Batch

%% MATLAB linear system creation using linearized PST results
tL = (0:0.01:15); % time to match PST d file time
modSig=zeros(1,size(tL,2)); % create blank mod signal same length as tL vector
modSig(find(tL> 1.0 ))= -0.025; % mirror logic from exciterModSig into input vector
modSig(find(tL>5))= 0; % mirror logic from exciterModSig into input vector

bsys = b_pm;
csys = [c_v;c_spd;c_pm];
G = ss(a_mat,bsys,csys,zeros(size(csys,1),size(bsys,2))); % create system using pst matricies

y = lsim(G,[modSig; zeros(1,length(modSig))],tL); % run input into state space system

% collect bus voltage magnitudes and adjust by initial conditions
linV = y(:,1:4)'; % rotate into col vectors

load PSTpath
load pMechStepNonLin

for busN = 1:size(linV,1)
    if strcmp('PSTv4', pstVer) || strcmp('pstSETO', pstVer)
        linV(busN,:) = linV(busN,:) + g.bus.bus(busN,2);
    else
        linV(busN,:) = linV(busN,:) + bus(busN,2);
    end
end

% collect machine speeds and adjust by initial condition
linSpd = y(:,5:6)'+ 1.0; % rotate into col vectors

% collect pm...
linPm = y(:,7:8)';% rotate to vector
if strcmp('PSTv4', pstVer) || strcmp('pstSETO', pstVer)
    linPm(1,:)= linPm(1,:)+ g.mac.pmech(1,1);
    linPm(2,:)= linPm(2,:)+ g.mac.pmech(2,1);
else
    linPm(1,:)= linPm(1,:) + pmech(1,1);
    linPm(2,:)= linPm(2,:) + pmech(2,1);
end

save linResults.mat tL linV linSpd modSig linPm

%% Clean up modulation file alterations.
% turbine governor pref mod
load PSTpath.mat
copyfile([PSTpath 'mpm_sig_ORIG.m'],[PSTpath 'mpm_sig.m']); % Replace original file

% Return new sub transient model back to OG
copyfile([PSTpath 'mac_sub_ORIG.m'],[PSTpath 'mac_sub.m']); % copy simulation specific data file to batch run location

%% temp file clean up
delete('PSTpath.mat')

%% load data for plot comparisons
load pMechStepNonLin.mat
load linResults.mat

%% compare Inputs to machine
figure
hold on
plot(tL,modSig)
if strcmp('PSTv4', pstVer) || strcmp('pstSETO', pstVer)
    plot(g.sys.t,g.mac.pm_sig(1,:),'--')
else
    plot(t,pm_sig(1,:),'--')
end
legend('Linear','Non-Linear','location','best')
title('Modulation Signal')

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
xlabel('Time [sec]')
ylabel('Voltage [PU]')

%% compare machine speeds
figure
hold on
legNames={};
for busN=1:size(linSpd,1)
    plot(tL,linSpd(busN,:))
    legNames{end+1}= ['Gen Speed ', int2str(busN), ' Linear'];
    
    if strcmp('PSTv4', pstVer) || strcmp('pstSETO', pstVer)
        plot(g.sys.t,g.mac.mac_spd(busN,:),'--')
    else
        plot(t,mac_spd(busN,:),'--')
    end
    
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
    
    if strcmp('PSTv4', pstVer) || strcmp('pstSETO', pstVer)
        plot(g.sys.t,g.mac.pmech(busN,:),'--')
    else
        plot(t,pmech(busN,:),'--')
    end
    
    legNames{end+1}= ['Gen Pm ', int2str(busN), ' non-Linear'];
end
legend(legNames,'location','best')
title('Machine Mechanical Power')
xlabel('Time [sec]')
ylabel('Mechanical Power [PU MW]')
