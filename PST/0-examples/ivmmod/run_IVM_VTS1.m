% Example of how to use ivmmod_dyn.m to model an IVM generator.
% Data file = d2m_ivmmod1.m
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
addpath([PSTpath, 'test', filesep]) % to handle new functionized code
save PSTpath.mat PSTpath
clear folderDepth pathParts pNdx PSTpath

%% Run nonlinear simulation and store results
clear all; close all; clc
load PSTpath
delete([PSTpath 'mac_sub.m']); 
copyfile([PSTpath 'mac_sub_NEW.m'],[PSTpath 'mac_sub.m']); %generator model

copyfile([PSTpath 'pss2.m'],[PSTpath 'pss.m']); % pss 2 model

delete([PSTpath 'DataFile.m']); 
copyfile('d2m_ivmmod1_VTS.m',[PSTpath 'DataFile.m']); %System data file

delete([PSTpath 'ivmmod_dyn.m']); 
copyfile('ivmmod_dyn_VTS1.m',[PSTpath 'ivmmod_dyn.m']); %Modulation file

% s_simu_Batch %Run PST
s_simu_BatchTestF %Run PST
% s_simu_BatchVTS %Run PST

save('ivmF'); %Save

%% restore default models
load PSTpath.mat
copyfile([PSTpath 'mac_sub_ORIG.m'],[PSTpath 'mac_sub.m']); % subtransient machine model
copyfile([PSTpath 'ivmmod_dyn_ORIG.m'],[PSTpath 'ivmmod_dyn.m']); %Modulation file
copyfile([PSTpath 'pss3.m'],[PSTpath 'pss.m']); % pss  model
delete PSTpath.mat

%% Plot
%load('Example1_NonlinearSim'); 

figure
subplot(411)
nb = 2; %Bus to plot
ng = 3; %Generator to plot (IVM generator at bus nb)
ni = 1; %IVM device (corresponds to generator ng and bus nb)

plot(g.sys.t,abs(g.bus.bus_v(nb,:)),'g');
hold on
plot(g.sys.t,g.ivm.ivmmod_e_sig(ni,:),'m','LineWidth',1.5);
plot(g.sys.t,g.mac.edprime(ng,:),'k--','LineWidth',1);
%ylim([1.1 1.25])
title('IVM Voltage')
legend(['bus ' num2str(nb)],'E_c','E','Location','SouthEast')
grid on
xlabel('Time [sec]')
ylabel('Voltage (pu)')

subplot(412)
plot(g.sys.t,g.ivm.ivmmod_d_sig(ni,:),'m','LineWidth',1.5);
hold on
plot(g.sys.t,g.mac.mac_ang(ng,:),'k--','LineWidth',1);
legend('{\delta}_c','{\delta}','Location','SouthEast')
grid on
xlabel('Time [sec]')
title('IVM Angle')
%ylim([0.1 0.35])
ylabel('IVM angle (rad)')

subplot(413)
plot(g.sys.t,g.mac.pelect(ng,:),'k','LineWidth',1.5)
grid on
xlabel('Time [sec]')
ylabel('P [MW PU]')
title('IVM Generated Real Power')

subplot(414)
plot(g.sys.t,g.mac.qelect(ng,:),'k','LineWidth',1.5)
grid on
xlabel('Time [sec]')
ylabel('Q [MVAR PU]')
title('IVM Generated Reactive Power')

%set(gcf,'Position',[520   0.5*378   560   1.5*420])

%% machine angle plots
figure
legNames={};
for n=1:size(g.mac.mac_ang,1)
    subplot(2,1,2)
    plot(g.sys.t, g.mac.mac_ang(n,:))
    hold on
    subplot(2,1,1)
    plot(g.sys.t, g.mac.mac_ang(n,:)-g.mac.mac_ang(1,:))
    hold on
    legNames{end+1}=['Machine ', int2str(n)];
end
subplot(2,1,2)
legend(legNames,'location','best')
grid on
ylabel('Angle [rad]')
xlabel('Time [sec]')
title('Machine Angles')

subplot(2,1,1)
legend(legNames,'location','best')
grid on
ylabel('Angle [rad]')
xlabel('Time [sec]')
title('Slack Referenced Machine Angles')