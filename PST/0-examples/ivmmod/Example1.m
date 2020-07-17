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
save PSTpath.mat PSTpath
clear folderDepth pathParts pNdx PSTpath

%% Run nonlinear simulation and store results
clear all; close all; clc
load PSTpath
delete([PSTpath 'mac_sub.m']); 
copyfile([PSTpath 'mac_sub_NEW.m'],[PSTpath 'mac_sub.m']); %generator model
delete([PSTpath 'DataFile.m']); 
copyfile('d2m_ivmmod1.m',[PSTpath 'DataFile.m']); %System data file
delete([PSTpath 'ivmmod_dyn.m']); 
copyfile('ivmmod_dyn_Example1.m',[PSTpath 'ivmmod_dyn.m']); %Modulation file

s_simu_Batch %Run PST
save('Example1_NonlinearSim'); %Save

%% Plot
%load('Example1_NonlinearSim'); 
figure
subplot(311)
nb = 2; %Bus to plot
ng = 3; %Generator to plot (IVM generator at bus nb)
ni = 1; %IVM device (corresponds to generator ng and bus nb)
plot(g.sys.t,abs(g.bus.bus_v(nb,:)),'k',g.sys.t,ivmmod_e_sig(ni,:),'r',g.sys.t,g.mac.edprime(ng,:),'b','LineWidth',2);
ylim([1.1 1.25])
legend(['bus ' num2str(nb)],'E_c','E','Location','SouthEast')
ylabel('Voltage (pu)')

subplot(312)
plot(g.sys.t,ivmmod_d_sig(ni,:),'r',g.sys.t,g.mac.mac_ang(ng,:),'b','LineWidth',2);
legend('{\delta}_c','{\delta}','Location','SouthEast')
ylim([0.1 0.35])
ylabel('IVM angle (rad)')

subplot(313)
plot(g.sys.t,g.mac.pelect(ng,:),'k','LineWidth',2)
ylabel('IVM real-power (pu)')

%set(gcf,'Position',[520   0.5*378   560   1.5*420])

%% clean up file manipulations
load PSTpath.mat
copyfile([PSTpath 'mac_sub_ORIG.m'],[PSTpath 'mac_sub.m']); % subtransient machine model
copyfile([PSTpath 'ivmmod_dyn_ORIG.m'],[PSTpath 'ivmmod_dyn.m']); %Modulation file
delete PSTpath.mat
