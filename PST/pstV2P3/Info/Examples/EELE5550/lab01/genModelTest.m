% Example of different generator models reaction to trip
% uses custom plotting function

clear all; close all; clc

%% Add pst path to MATLAB
% generate relative path generically
folderDepth = 4; % depth of current directory from main PST directory

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
delete([PSTpath 'DataFile.m']); copyfile('dlab_01_classical.m',[PSTpath 'DataFile.m']); %System data file
delete([PSTpath 'pwrmod_dyn.m']); %copyfile('pwrmod_dyn_Example2.m',[PSTpath 'pwrmod_dyn.m']); %Modulation file
s_simu_Batch %Run PST <- this is the main file to look at for simulation workings
save('gen1_classical','t','bus_v','pelect','mac_spd'); %Save t and bus_v results

%% Run nonlinear simulation and store results
clear all; close all; clc
load PSTpath
delete([PSTpath 'DataFile.m']); copyfile('dlab_01_2axis.m',[PSTpath 'DataFile.m']); %System data file
delete([PSTpath 'pwrmod_dyn.m']); %copyfile('pwrmod_dyn_Example2.m',[PSTpath 'pwrmod_dyn.m']); %Modulation file
s_simu_Batch %Run PST <- this is the main file to look at for simulation workings
save('gen1_2axis','t','bus_v','pelect','mac_spd'); %Save t and bus_v results

%% Run nonlinear simulation and store results
clear all; close all; clc
load PSTpath
delete([PSTpath 'DataFile.m']); copyfile('dlab_01_subt.m',[PSTpath 'DataFile.m']); %System data file
delete([PSTpath 'pwrmod_dyn.m']); %copyfile('pwrmod_dyn_Example2.m',[PSTpath 'pwrmod_dyn.m']); %Modulation file
s_simu_Batch %Run PST <- this is the main file to look at for simulation workings
save('gen1_subt','t','bus_v','pelect','mac_spd'); %Save t and bus_v results

%% Plotting
title_str = 'Comparison of Generator Models';
plot_end = 9;

stability_plot( 'gen1_classical.mat', 'Classical', ...
    'gen1_2axis.mat', 'Two-Axis', ...
    'gen1_subt.mat', 'Sub-Transient', ...
    title_str, plot_end )

%% temp file clean up
delete('PSTpath.mat')
delete('sim_fle.mat')