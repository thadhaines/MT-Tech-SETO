% Example of how to use ivmmod_dyn.m to simulate voltage source behind impedance
% thru a model

% Reused code from pwrmod used to test ivmod
% Data file = d2m_ivmmod1.m - linear method using ivmod
% tested as working in octave 5/12/20

% CONTAINS LEFT OVER CODE THAT IS COMMENTED OUT 
% WIP example based on pwrmod
% Data file = d2m_pwrmod2.m - non-linear method from previous example

clear all; close all; clc

%% Add pst path to MATLAB
% generate relative path generically
folderDepth = 3; % depth of current directory from main PST directory

pathParts = strsplit(pwd, filesep);
PSTpath = pathParts(1);

for pNdx = 2:max(size(pathParts))-folderDepth
    PSTpath = [char(PSTpath), filesep, char(pathParts(pNdx))];
end
PSTpath = [char(PSTpath), filesep];

addpath(PSTpath)
save PSTpath PSTpath

%% Run nonlinear simulation and store results
clear all; %close all; clc
load PSTpath
delete([PSTpath 'DataFile.m']); copyfile('d2m_ivmmod1.m',[PSTpath 'DataFile.m']); %System data file
delete([PSTpath 'ivmmod_dyn.m']); copyfile('ivmmod_dyn_Example1.m',[PSTpath 'ivmmod_dyn.m']); %Modulation file
s_simu_Batch %Run PST
save('Exmaple1_NonlinearSim.mat'); %Save

%% clean up
delete('PSTpath') % generated in this script
delete('PSTpath.mat') % generated in this script
delete('sim_fle') % generated from PST...
delete('sim_fle.mat') % generated from PST...

%% Plot results from earlier sim
% load trip on bus 6, gen 2 = ivm, gen1 hydrogov , 
figure
plot(t, ivmmod_e_sig)
hold on
plot(t, abs(bus_v(2,:)))
legend('Input Signal','Bus Voltage')
title('ivmmod e sig')
figure
plot(t, ivmmod_d_sig)
hold on
plot(t, angle(bus_v(2,:)))
legend('Input Signal','Bus Voltage Angle')
title('ivmmod d sig')


%{
 Left over/ reused code from pwrmod

%% Build Linear model, simulate, and store results
%Build PST linear model with real-power input and voltage mag and angle
%outputs
clear all; %close all; clc
load PSTpath
delete([PSTpath 'DataFile.m']); 
copyfile('d2m_pwrmod2.m',[PSTpath 'DataFile.m']); %System data file
delete([PSTpath 'pwrmod_dyn.m']); 
copyfile('pwrmod_dyn_Example2.m',[PSTpath 'pwrmod_dyn.m']); %Modulation file


svm_mgen_Batch %Conduct linearization
Gpst = ss(a_mat,b_pwrmod_p,[c_v([2 3],:); c_ang([2 3],:)],zeros(4,2));

%Build pwrmod_dyn linear model 
%Inputs are Pref1, abs(bus_v(2,:)), Pref2, abs(bus_v(3,:))
%Ouputs are Ipcmd at bus 2, Ipcmd at bus 3
Tpord = [0.25; 0.35]; %Power order filter time constant
Tv = [0.05; 0.15]; %Voltage measurment filter
A = zeros(4,4); A(1,1) = -1/Tpord(1); A(2,2) = -1/Tv(1); A(3,3) = -1/Tpord(2); A(4,4) = -1/Tv(2); 
B = eye(4,4); B(1,1) = 1/Tpord(1); B(2,2) = 1/Tv(1); B(3,3) = 1/Tpord(2); B(4,4) = 1/Tv(2); 
C = zeros(2,4); C(1,1) = 1/bus(2,2); C(2,3) = 1/bus(3,2); %Eg, linearize Ipcmd(1) = x(1)/x(2)
D = zeros(2,4);
Gsolar = ss(A,B,C,D);
clear A B C D

%Connect the two linear models together
GopenLoop = series(Gsolar,Gpst);
G = feedback(GopenLoop,eye(2,2),[2;4],[1;2],1);
clear Gsolar Gpst GopenLoop

%Simulate linear model
tL = [0:1/120:10]';
N = length(tL);
Pref = [-0.0001*(stepfun(tL,1)-stepfun(tL,1.5)) zeros(N,1) 0.0002*(stepfun(tL,4)-stepfun(tL,4.5)) zeros(N,1) ]; %Pref pulses
load('Exmaple2_NonlinearSim','t','bus_v');
y = lsim(G,Pref,tL);
v = y(:,1:2) + ones(N,1)*abs(bus_v([2 3],1))';
a = y(:,3:4) + ones(N,1)*angle(bus_v([2 3],1))';
bus_vL = transpose(v.*exp(1i*a));
tL = tL';
save('Exmaple2_LinearSim','tL','bus_vL'); %Save linear results

%% Plot linear vs nonlinear
clear all; close all; clc
load('Exmaple2_NonlinearSim','t','bus_v','pwrmod_p_st','pwrmod_q_st');
load('Exmaple2_LinearSim','tL','bus_vL');
figure
subplot(411)
nb = 2; %Bus to plot
plot(t,abs(bus_v(nb,:)),'k',tL,abs(bus_vL(1,:)),'r');
ylabel(['bus ' num2str(nb) ' V (abs)'])
f = 1e3*angle(bus_v(nb,2:end)./bus_v(nb,1:end-1))./(2*pi*diff(t)); f = [f f(end)];
fL = 1e3*angle(bus_vL(1,2:end)./bus_vL(1,1:end-1))./(2*pi*diff(tL)); fL = [fL fL(end)];
subplot(412)
plot(t,f,'k',tL,fL,'r')
ylabel(['bus ' num2str(nb) ' (mHz)'])
subplot(413)
nb = 3;
plot(t,abs(bus_v(nb,:)),'k',tL,abs(bus_vL(2,:)),'r');
ylabel(['bus ' num2str(nb) ' V (abs)'])
f = 1e3*angle(bus_v(nb,2:end)./bus_v(nb,1:end-1))./(2*pi*diff(t)); f = [f f(end)];
fL = 1e3*angle(bus_vL(2,2:end)./bus_vL(2,1:end-1))./(2*pi*diff(tL)); fL = [fL fL(end)];
subplot(414)
plot(t,f,'k',tL,fL,'r')
ylabel(['bus ' num2str(nb) ' (mHz)'])
set(gcf,'Position',[360 202 560 720]);

%% Plot nonlinear P and Q injected power
figure
subplot(411)
plot(t,pwrmod_p_st(1,:),'k')
ylabel('Bus 2 I_P (pu)')
subplot(412)
plot(t,pwrmod_p_st(2,:),'k')
ylabel('Bus 3 I_P (pu)')
subplot(413)
plot(t,pwrmod_q_st(1,:),'k')
ylabel('Bus 2 I_Q (pu)')
subplot(414)
plot(t,pwrmod_q_st(2,:),'k')
ylabel('Bus 3 I_Q (pu)')
xlabel('Time (sec.)')
set(gcf,'Position',[360 202 560 720]);
%}
