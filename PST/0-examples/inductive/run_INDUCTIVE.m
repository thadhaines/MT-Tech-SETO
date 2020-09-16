% inductive machine/generator and inductive load case
% Tested as working in all Versions
% output the same in all versions (iff all models are the same)

clear all; close all; clc
case2run = 2; % 1= fault, 2 = load step

%% Add pst path to MATLAB
% generate relative path generically
folderDepth = 2; % depth of current directory from main PST directory
pstVer =   'PSTv4'; % 'pstSETO'; %  'pstV2p3';%  'pstV3P1';% 
pathParts = strsplit(pwd, filesep);
PSTpath = pathParts(1);

for pNdx = 2:max(size(pathParts))-folderDepth
    PSTpath = [char(PSTpath), filesep, char(pathParts(pNdx))];
end
PSTpath = [char(PSTpath), filesep, pstVer, filesep];

addpath(PSTpath)

save PSTpath.mat PSTpath pstVer case2run
clear folderDepth pathParts pNdx PSTpath

%% Run nonlinear simulation and store results
clear all; clc
load PSTpath.mat

if case2run == 1
    copyfile('data3Trip.m',[PSTpath 'DataFile.m']); % trip case
else
    copyfile('data3mIg.m',[PSTpath 'DataFile.m']); % load step case
    
    if strcmp('PSTv4', pstVer) || strcmp('pstSETO', pstVer)
        copyfile( 'ml_sig_smallStepG.m',[PSTpath 'ml_sig.m']); % For global g
    else
        copyfile( 'ml_sig_smallStep.m',[PSTpath 'ml_sig.m']); % for v 2.3 and 3.1
    end
    
end

copyfile([PSTpath 'mac_ind2.m'],[PSTpath 'mac_ind.m']); % ensure newest model

livePlotFlag = 1;
if strcmp('PSTv4', pstVer)
    copyfile([PSTpath 'livePlot_1.m'],[PSTpath 'livePlot.m']); % specify plot operation
    s_simu
    copyfile([PSTpath 'livePlot_ORIG.m'],[PSTpath 'livePlot.m']); % reset live plot
else
    s_simu_Batch %Run PST non-linear sim
end

% reset modulation file
copyfile([PSTpath 'ml_sig_ORIG.m'],[PSTpath 'ml_sig.m']); % copy system data file to batch run location

%% Save cleaned output data
save([pstVer,'INDnonLIN.mat']); %Save simulation outputs

if case2run == 2 % only execute linear sim if load step
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
    load PSTpath
    for busN = 1:size(linV,1)
        if strcmp('PSTv4', pstVer) || strcmp('pstSETO', pstVer)
            linV(busN,:) = linV(busN,:) + g.bus.bus(busN,2);
            linAng(busN,:) = linAng(busN,:) + deg2rad(g.bus.bus(busN,3));
        else
            linV(busN,:) = linV(busN,:) + bus(busN,2);
            linAng(busN,:) = linAng(busN,:) + deg2rad(bus(busN,3));
        end
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
    if strcmp('PSTv4', pstVer) || strcmp('pstSETO', pstVer)
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
    
    %% compare bus angle
    figure
    hold on
    legNames={};
    for busN=1:size(linAng,1)
        plot(tL,linAng(busN,:))
        legNames{end+1}= ['Bus ', int2str(busN), ' Linear'];
        if strcmp('PSTv4', pstVer) || strcmp('pstSETO', pstVer)
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
end % end linear comparison
