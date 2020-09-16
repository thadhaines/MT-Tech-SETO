% Example of two machine system with exciter modulation
% smaller system used to more easily study inner workings of PST
% One steam gov gen, One hydro gov.
% only 1 exciter
% Hydro gov twice size as Steam.

% NOTE: live plotting disabled in d_ file to test speed up and better compare MATLAB to Octave performance
% (Octave seems a bit faster when plotting is removed from simulation run)

%% Run nonlinear simulation and store results
clear all; close all; clc
load PSTpath.mat
delete([PSTpath 'DataFile.m']); % ensure batch datafile is cleared
copyfile('d_EXC_type0.m',[PSTpath 'DataFile.m']); % copy system data file to batch run location

% Handle exciter modulation file placement etc.
%copyfile([PSTpath 'mexc_sig.m'],[PSTpath 'mexc_sig_ORIG_TMP.m']); % save copy of original ml_sig file
delete([PSTpath 'mexc_sig.m']); % ensure ml_sig file is empty
copyfile('exciterModSigG.m',[PSTpath 'mexc_sig.m']); % copy simulation specific data file to batch run location

s_simu %Run PST non-linear sim

%% Simulation variable cleanup
% Clear any varables that contain only zeros
varNames = who()'; % all variable names in workspace
clearedVars = {}; % cell to hold names of deleted 'all zero' variables

for vName = varNames
    try
        zeroTest = eval(sprintf('all(%s(:)==0)', vName{1})); % check if all zeros
        if zeroTest
            eval(sprintf('clear %s',vName{1}) ); % clear variable
            clearedVars{end+1} = vName{1}; % add name to cell for reference
        end
    catch ME
        disp(ME.message)
        disp(vName)
    end
    
end
clear varNames vName zeroTest

%% Gernalize modulation input signal storage
modSigNL = g.exc.exc_sig;

%% Save cleaned output data
save('excT0NL.mat'); %Save all non-linear simulation outputs

%% PST linear system creation
clear all; close all;
svm_mgen_Batch

%% MATLAB linear system creation using linearized PST results
tL = (0:0.01:5); % time to match PST d file time
modSig=zeros(1,size(tL,2)); % create blank mod signal same length as tL vector
modSig(find(tL>0.5))= 0.01; % mirror logic from exciterModSig into input vector

bsys = b_vr; % input to system is exciter voltage reference
csys = [c_v;c_spd;c_pm];
G = ss(a_mat,bsys,csys,zeros(size(csys,1),size(bsys,2))); % create system using pst matricies

y = lsim(G,modSig,tL); % run input into state space system

% collect bus voltage magnitudes and adjust by initial conditions
linV = y(:,1:4)'; % rotate into col vectors
for busN = 1:size(linV,1)
    linV(busN,:) = linV(busN,:) + g.bus.bus(busN,2);
end

% collect machine speeds and adjust by initial condition
linSpd = y(:,5:6)'+1.0; % rotate into col vectors

% collect pm.
linPm = y(:,7:8)';% rotate to vector
linPm(1,:)= linPm(1,:)+ g.mac.pmech(1,1);
linPm(2,:)= linPm(2,:)+ g.mac.pmech(2,1);
save excT0LIN.mat tL linV linSpd linPm modSig
clear all

%% Clean up modulation file alterations...
load PSTpath.mat
delete([PSTpath 'mexc_sig.m']); % remove simulation specific ml_sig file
copyfile([PSTpath 'mexc_sig_ORIG.m'],[PSTpath 'mexc_sig.m']); % Replace original file
