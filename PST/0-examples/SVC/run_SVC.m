% svc test
% Tested as working in all Versions
% output the same in all versions
clear all; close all; clc

%% Add pst path to MATLAB
% generate relative path generically
folderDepth = 2; % depth of current directory from main PST directory
pstVer =  'PSTv4'; % 'pstSETO'; %  'pstV2p3';%  'pstV3P1';% 
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
copyfile('d2a_dceREF.m',[PSTpath 'DataFile.m']); % copy system data file to batch run location

% move modulation file
if strcmp('PSTv4', pstVer) || strcmp('pstSETO', pstVer)    
    copyfile('msvc_sig_SmallStepG.m',[PSTpath 'msvc_sig.m']);     
else
    copyfile('msvc_sig_SmallStep.m',[PSTpath 'msvc_sig.m']); 
end

if strcmp('PSTv4', pstVer)
    s_simu
else
    s_simu_Batch %Run PST
end

% reset modulation file
copyfile([PSTpath 'msvc_sig_ORIG.m'],[PSTpath 'msvc_sig.m']); % copy system data file to batch run location

%% Save output data
save([pstVer,'SVCnonLIN.mat']); %Save simulation outputs

%% PST linear system creation
clear all; close all;

pssGainFix = 1;
svm_mgen_Batch

%% MATLAB linear system creation using linearized PST results
tL = (0:0.001:10);              % time to match PST d file time
modSig=zeros(1,size(tL,2));     % create blank mod signal same length as tL vector
modSig(find(tL>= 0.5 ))= 0.01;  % mirror logic from exciterModSig into input vector
modSig(find(tL>= 1.5))= 0;      % mirror logic from exciterModSig into input vector

bsys = b_svc;
csys = [c_v;c_spd;c_pm];
G = ss(a_mat,bsys,csys,zeros(size(csys,1),size(bsys,2))); % create system using pst matricies

y = lsim(G,modSig,tL); % run input into state space system

% collect bus voltage magnitudes and adjust by initial conditions
linV = y(:,1:size(c_v,1))'; % rotate into col vectors
load PSTpath
for busN = 1:size(linV,1)
    if strcmp('PSTv4', pstVer) || strcmp('pstSETO', pstVer)
        linV(busN,:) = linV(busN,:) + g.bus.bus(busN,2);
    else
        linV(busN,:) = linV(busN,:) + bus(busN,2);
    end
end

% collect machine speeds and adjust by initial condition
spdStart = size(c_v,1)+1;
spdEnd = spdStart + size(c_spd,1)-1;
linSpd = y(:,spdStart:spdEnd )'+ 1.0; % rotate into col vectors

% collect pm...
pmStart = spdEnd+1;
linPm = y(:,pmStart:end)';% rotate to vector

% required adjustments
name = [pstVer,'SVCnonLIN.mat'];
feval('load', name)

for pmAdj = 1:size(linPm,1)
    if strcmp('PSTv4', pstVer) || strcmp('pstSETO', pstVer)
        linPm(pmAdj,:)= linPm(pmAdj,:)+ g.mac.pmech(pmAdj,1);
    else
        linPm(pmAdj,:)= linPm(pmAdj,:)+ pmech(pmAdj,1);
    end
end

save linResults.mat tL linV linSpd modSig linPm

%% plot comparisons
name = [pstVer,'SVCnonLIN.mat'];
feval('load', name)
load linResults.mat

%% compare mod inputs
figure
hold on
plot(tL,modSig)
if strcmp('PSTv4', pstVer) || strcmp('pstSETO', pstVer)
    plot(g.sys.t,g.svc.svc_sig,'--')
else
    plot(t,svc_sig,'--')
end
legend('Linear','Non-Linear','location','best')
title('Governor Pref Modulation Signal')


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

%% temp file clean up
delete('PSTpath.mat')