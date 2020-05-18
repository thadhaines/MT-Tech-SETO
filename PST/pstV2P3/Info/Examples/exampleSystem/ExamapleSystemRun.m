% Example of 2 cycle 3 phase line fault using system defined in PSTV3
% manual.
% Many Plots of PST data provided at end

% s_simu_Batch modified to show warnings at various simulation points - the
% warning also supplies line numbers which are useful

% when k == 50, i.e. the 50th simulation loop interation, any call to the
% network solution i_simu is identified with a warning. This shows that the
% network is solved twice... 

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
clear all; close all; clc
load PSTpath
delete([PSTpath 'DataFile.m']); % ensure batch datafile is cleared
copyfile('d_exampleSystem3pFault.m',[PSTpath 'DataFile.m']); %System data file
s_simu_Batch %Run PST <- this is the main file to look at for simulation workings
save('exampleSystem01'); %Save t and bus_v results

%% temp file clean up
delete('PSTpath.mat')
delete('sim_fle.mat')

%% Plotting init
clear all % ensure only saved data is plotted
load('exampleSystem01.mat')

%% Plotting of PST outputs
% Synchronous generator info
figure
plot(t,pmech)
title([{'pmech'}; {'generator mechanical power output'}])
xlabel('Time [sec]')
ylabel('MW [PU]')

figure
plot(t,pelect)
title([{'pelect'}; {'generator active power output'}])
xlabel('Time [sec]')
ylabel('MW [PU]')

figure
plot(t,qelect)
title([{'qelect'}; {'generator reactive power output'}])
xlabel('Time [sec]')
ylabel('MVAR [PU]')

figure
plot(t,mac_spd)
title([{'mac\_spd'}; {'machine speed'}])
xlabel('Time [sec]')
ylabel('Speed [PU]')

figure
plot(t,mac_ang)
title([{'mac\_ang'}; {'machine angle'}])
xlabel('Time [sec]')
ylabel('Angle [rads]')

% zeros: pm_sig, mac_ref, 
% used files: tg 1-3, pss 1-3 (states)

% Exciter info
figure
plot(t,Efd)
title([{'Efd'}; {'excitation output voltage (field voltage)'}])
xlabel('Time [sec]')
ylabel('Voltage [PU]')

% bus info
figure
plot(t,abs(bus_v))
title([{'abs(bus\_v)'}; {'bus voltage magnitude'}])
xlabel('Time [sec]')
ylabel('Voltage [PU]')

figure
plot(t,angle(bus_v))
title([{'angle(bus\_v)'}; {'bus voltage angle'}])
xlabel('Time [sec]')
ylabel('Angle [rads]')

% % Branch current - Probably just current phasor
% figure
% plot(t,real(ilf))
% title([{'real(ilf)'}; {'active current line flow'}])
% xlabel('Time [sec]')
% ylabel('Current [PU]')
% 
% figure
% plot(t,imag(ilf))
% title([{'imag(ilf)'}; {'reactive current line flow'}])
% xlabel('Time [sec]')
% ylabel('Current [PU]')

%% end comments
%{
Current obvious issues:
exclusive use of global variables
switching of fault status via replacing Y matricies requires further study
batch run assumes 60 Hz, 100 MVA base
no identifier for multiple lines connecting to same bus
no string identifier for any object
no running log of bus load powers
seperation of data from object identifier -> labeling/identification requires cross referencing
matrix definitions get large/clunky/unweildy
bus_v written multiple times per time step (shown if DEBUG = 1 in d_ file)
i_simu called multiple times per time step
time step larger than 1 cycle causes system to 'blow up' i.e. too large of an integration step
anlges will require a reference to the slack bus angle and be adjusted appropriately to produce 'standard reference' plot
%}