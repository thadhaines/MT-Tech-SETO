% Example of how to use pwrmod_dyn.m to modulate the power into a bus.
% Data file = d2m_pwrmod1.m

%% Octave specific
% https://wiki.octave.org/Differences_between_Octave_and_Matlab
warning('off', 'Octave:possible-matlab-short-circuit-operator'); # supress warning about | and & usage vs || and &&
pkg load control % for ss functionality

% Tested with octave 5.2.0
% requires stepfun for linear model - this function is recreated from matlab function
% file saving / loading requires explicit .mat specification - sometimes?

%% script start
close all; clear; clc

%% Add pst path 
% generate relative path generically
folderDepth = 3; % depth of current directory from main PST directory

pathParts = strsplit(pwd, filesep);
PSTpath = pathParts(1);

for pNdx = 2:max(size(pathParts))-folderDepth
    PSTpath = [char(PSTpath), filesep, char(pathParts{pNdx})];
end
PSTpath = [char(PSTpath), filesep];

addpath(PSTpath)
save PSTpath PSTpath

%% Run nonlinear simulation and store results
clear all; close all; clc
load PSTpath
delete([PSTpath 'DataFile.m']); % clear previous system settings (if applicable)
copyfile('d2m_pwrmod1.m',[PSTpath 'DataFile.m']); % Place new system data file
delete([PSTpath 'pwrmod_dyn.m']); % clear previous pwrmod settings (if applicable)
copyfile('pwrmod_dyn_Example1.m',[PSTpath 'pwrmod_dyn.m']); % Place new modulation file
s_simu_Batch %Run PST <- this is the main file to look at for simulation workings
save('Exmaple1_NonlinearSim.mat','t','bus_v','pwrmod_p_st','pwrmod_q_st'); %Save t and bus_v results

%% Build linear model, simulate, and store results
%Build linear model
clear all; clc;
load PSTpath
delete([PSTpath 'DataFile.m']); % clear previous system settings (if applicable)
copyfile('d2m_pwrmod1.m',[PSTpath 'DataFile.m']); %System data file
delete([PSTpath 'pwrmod_dyn.m']); % clear previous pwrmod settings (if applicable)
copyfile('pwrmod_dyn_Example1.m',[PSTpath 'pwrmod_dyn.m']); %Modulation file
svm_mgen_Batch %Conduct linearization
save('Example1_Linear.mat','a_mat','b_pwrmod_p','c_v','c_ang');

%% Plot nonlinear P and Q injected power
load('Exmaple1_NonlinearSim.mat')
load PSTpath
figure(1)
subplot(411)
plot(t,pwrmod_p_st(1,:),'k')
ylabel('Bus 2 P (pu)')
title('Real Power Injection')
subplot(412)
plot(t,pwrmod_p_st(2,:),'k')
ylabel('Bus 3 P (pu)')
subplot(413)
plot(t,pwrmod_q_st(1,:),'k')
ylabel('Bus 2 Q (pu)')
subplot(414)
plot(t,pwrmod_q_st(2,:),'k')
ylabel('Bus 3 Q (pu)')
xlabel('Time (sec.)')
%set(gcf,'Position',[360 202 560 720]);

%% if using statespace from Graham: shadows built in functions... -> breaks 
% i.e. don't do this
% ssPath = [char(PSTpath), 'statespace', filesep, '@stsp', filesep]
% addpath(ssPath)

%% Simulate linear model
Gv = ss(a_mat,b_pwrmod_p,c_v,zeros(6,2));
Ga = ss(a_mat,b_pwrmod_p,c_ang,zeros(6,2));
tL = [0:1/120:10]';
u = [-0.0001*(stepfun(tL,1)-stepfun(tL,1.5)) 0.0002*(stepfun(tL,4)-stepfun(tL,4.5))]; 
load('Exmaple1_NonlinearSim','t','bus_v');
v = lsim(Gv,u,tL);
v = v + ones(size(v,1),1)*abs(bus_v(1:6,1))';
a = lsim(Ga,u,tL);
a = a + ones(size(a,1),1)*angle(bus_v(1:6,1))';
bus_vL = transpose(v.*exp(1i*a));
tL = tL';
save('Exmaple1_LinearSim','tL','bus_vL'); %Save linear results

%% Plot Nonlinear vs Linear
%clear all; clc
load('Exmaple1_NonlinearSim','t','bus_v');
load('Exmaple1_LinearSim','tL','bus_vL');
figure(2)
subplot(411)
nb = 2; %Bus to plot
plot(t,abs(bus_v(nb,:)),'k',tL,abs(bus_vL(nb,:)),'r');
title('Voltage and Frequency Comparison - Power Injection')
ylabel(['bus ' num2str(nb) ' V (abs)'])
f = 1e3*angle(bus_v(nb,2:end)./bus_v(nb,1:end-1))./(2*pi*diff(t)); 
f = [f f(end)];
fL = 1e3*angle(bus_vL(nb,2:end)./bus_vL(nb,1:end-1))./(2*pi*diff(tL)); 
fL = [fL fL(end)];
subplot(412)
plot(t,f,'k',tL,fL,'r')
ylabel(['bus ' num2str(nb) ' (mHz)'])
subplot(413)
nb = 3; %Bus to plot
plot(t,abs(bus_v(nb,:)),'k',tL,abs(bus_vL(nb,:)),'r');
ylabel(['bus ' num2str(nb) ' V (abs)'])
f = 1e3*angle(bus_v(nb,2:end)./bus_v(nb,1:end-1))./(2*pi*diff(t)); 
f = [f f(end)];
fL = 1e3*angle(bus_vL(nb,2:end)./bus_vL(nb,1:end-1))./(2*pi*diff(tL)); 
fL = [fL fL(end)];

subplot(414)
plot(t,f,'k',tL,fL,'r')
ylabel(['bus ' num2str(nb) ' (mHz)'])
legend({'Non-Linear','Linear'},'location','best')

%set(gcf,'Position',[360 202 560 720]);