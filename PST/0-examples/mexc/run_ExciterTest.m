% Example of two machine system with exciter modulation
% smaller system used to more easily study inner workings of PST
% One steam gov gen, One hydro gov.
% only 1 exciter
% Hydro gov twice size as Steam.
% Works in all PST versions

clear all; close all; clc

%% Add pst path to MATLAB
% generate relative path generically
folderDepth = 2; % depth of current directory from main PST directory
pstVer =   'PSTv4'; %   'pstSETO'; %  'pstV2P3'; % 'pstV3p1'; %
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
copyfile('d_OMIB_EXC.m',[PSTpath 'DataFile.m']); % copy system data file to batch run location

% Handle exciter modulation file placement etc.
if strcmp('PSTv4', pstVer) || strcmp('pstSETO', pstVer)
    copyfile('exciterModSigG.m',[PSTpath 'mexc_sig.m']);
else
    copyfile('exciterModSig.m',[PSTpath 'mexc_sig.m']);
end

if strcmp('PSTv4', pstVer)
    s_simu
else
    s_simu_Batch %Run PST non-linear sim
end


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
load PSTpath
for busN = 1:size(linV,1)
    if strcmp('PSTv4', pstVer) || strcmp('pstSETO', pstVer)
        linV(busN,:) = linV(busN,:) + g.bus.bus(busN,2);
else
        linV(busN,:) = linV(busN,:) + bus(busN,2);
end

end

% collect machine speeds and adjust by initial condition
linSpd = y(:,5:6)'+1.0; % rotate into col vectors
save linResults.mat tL linV linSpd excSig
clear all

%% Clean up modulation file alterations...
load PSTpath.mat

copyfile([PSTpath 'mexc_sig_ORIG.m'],[PSTpath 'mexc_sig.m']); % Replace original file

%% temp file clean up
delete('PSTpath.mat')

%% plot comparisons
load exciterTest.mat
load linResults.mat
%% compare Inputs to exciter
figure
hold on
plot(tL,excSig)
    if strcmp('PSTv4', pstVer) || strcmp('pstSETO', pstVer)
        plot(g.sys.t,g.exc.exc_sig,'--')
else
        plot(t,exc_sig,'--')
    end



legend('Linear','Non-Linear','location','best')
title('Exciter Modulation Signal')
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
    plot(g.sys.t,g.mac.mac_spd(busN,:),'--')
else
    plot(t,mac_spd(busN,:),'--')
            end
        
    legNames{end+1}= ['Gen Speed ', int2str(busN), ' non-Linear'];
    
end
legend(legNames,'location','best')
title('Machine Speeds')