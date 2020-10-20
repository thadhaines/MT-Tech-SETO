% Run 3AreaTest
% Only functional in PSTv4

clear; close all; clc

%% Add correct PST verstion path to MATLAB in a generic way
folderDepth = 2; % depth of current directory from main PST directory
pstVer =   'PSTv4';
pathParts = strsplit(pwd, filesep);
PSTpath = pathParts(1);

for pNdx = 2:max(size(pathParts))-folderDepth
    PSTpath = [char(PSTpath), filesep, char(pathParts(pNdx))];
end
PSTpath = [char(PSTpath), filesep, pstVer, filesep];

addpath(PSTpath)
save PSTpath.mat PSTpath pstVer
clear folderDepth pathParts pNdx PSTpath scenario

%% Run nonlinear simulation and store results
clear all; close all; clc
load PSTpath.mat

copyfile('3AreaTest.m',[PSTpath 'DataFile.m']); % system data file
copyfile([PSTpath 'livePlot_2.m'],[PSTpath 'livePlot.m']); % set live plot
%copyfile('ml_sig_smallStepG.m',[PSTpath 'ml_sig.m']); % load step
copyfile('ml_sig_smallStep_fdepG.m',[PSTpath 'ml_sig.m']); % frequency dependent loads and load step

%livePlotFlag = 1; % not ideal for extended term sim - may cause crashes/freezing

s_simu

%% Clean up modulation file alterations.
copyfile([PSTpath 'ml_sig_ORIG.m'],[PSTpath 'ml_sig.m']); % load step
copyfile([PSTpath 'livePlot_ORIG.m'],[PSTpath 'livePlot.m']); % live plot

%% Save cleaned output data
save('xxx.mat'); %Save simulation outputs

%% temp file clean up
delete('PSTpath.mat')

%% ===============================================================

%% Plots
xevents = g.sys.sw_con(:,1)';
plotCell = { ...
    %f1, f2
    'mac','pmech';
   % 'tg','tg_sig';
    'mac','mac_spd';
    'mac','pelect';
%    'mac','qelect';
%     'mac','cur_re';
%     'mac','cur_im';
%     'mac','mac_ang';
%     'mac','psi_re';
%     'mac','psi_im';
    };

% nS = find(g.sys.t > 24);
% nS = nS(1);
% nE = find(g.sys.t > 28);
% nE = nE(1);

axLim = [24,28];
lnClr=[0,0,0; 0.66, 0.66, 0.66; 1, 0, 1]; % for custom colors

figure % to plot all in one plot
for n = 1:size(plotCell,1)
    f1 = plotCell{n,1};
    f2 = plotCell{n,2};
    %figure
    subplot(size(plotCell,1),1, n)
    
    for m = 1:size(g.(f1).(f2),1)
        if m<=3
            plot(g.sys.t, g.(f1).(f2)(m,:),'color',lnClr(m,:),'linewidth',1.25)
        else
            plot(g.sys.t, g.(f1).(f2)(m,:),'color','linewidth',1.25)
        end
        hold on
    end
    
    xlabel('Time [sec]')
    % handle multiple y labels
    if strcmp(f2, 'pmech') || strcmp(f2, 'pelect')
        ylabel('MW [PU]')
    elseif strcmp(f2, 'tg_sig')
        ylabel('Mod Signal [PU]')
    elseif strcmp(f2, 'mac_spd')
        ylabel('Speed [PU]')
    elseif strcmp(f2, 'qelect')
        ylabel('MVAR [PU]')
    end
    set(gca, 'XTick', xevents);
    legend({'Gen 1', 'Gen 2', 'Gen 3'},'location','best')
    
    %xlim(axLim) % for detail plots
    
    title([f1,'.',f2], 'Interpreter','None')
    grid on
end

%% AGC signals
figure
plot(g.sys.t, real(g.area.area(1).icA-g.area.area(1).icS)*g.sys.basmva,'linewidth',2)
hold on
plot(g.sys.t, real(g.area.area(2).icA-g.area.area(2).icS)*g.sys.basmva,'--','linewidth',2)
plot(g.sys.t, real(g.area.area(3).icA-g.area.area(3).icS)*g.sys.basmva,':','linewidth',2)
legend({'Area 1','Area 2', 'Area 3'}, 'location', 'best')
title('Interchange Error')
xlabel('Time [seconds]')
ylabel('IC Error [MW]')
grid on

%% frequency dependant load mod
figure
subplot(2, 1, 1)
plot(g.sys.t, g.lmod.lmod_sig(1,:))
legend('Area 1')
title('Load Modulation')
xlabel('Time [seconds]')
ylabel('Load Mod. [PU]')
grid on

subplot(2, 1, 2)
plot(g.sys.t, g.lmod.lmod_sig(2,:))
hold on
plot(g.sys.t, g.lmod.lmod_sig(3,:))
legend('Area 2','Area 3')
xlabel('Time [seconds]')
ylabel('Load Mod. [PU]')
grid on

%% System Bus voltages
% figure
% legNames = {};
% for n=1:size(g.bus.bus_v,1)-1
%     plot(g.sys.t, abs(g.bus.bus_v(n,:)))
%     legNames{end+1} = ['Bus ',num2str(g.bus.bus(n,1))];
%     hold on
% end
% 
% title('All Sustem Bus Voltages')
% ylabel('Bus Voltage [PU]')
% xlabel('Time [seconds]')
% set(gca, 'XTick', xevents);
% legend(legNames, 'location','best')
% grid on

%% Line currents feeding load
% figure
% legNames = {};
% for n=1:g.lmon.n_lmon
%     subplot(2,1,1)
%     plot(g.sys.t, real(g.lmon.line(n).sFrom),'color',lnClr(n,:),'linewidth',1.25)
%     hold on
%     grid on
%     
%     subplot(2,1,2)
%     plot(g.sys.t, imag(g.lmon.line(n).sFrom),'color',lnClr(n,:),'linewidth',1.25)
%     hold on
%     grid on
%     
%     legNames{end+1} = ['Line ',num2str(g.lmon.line(n).FromBus), ' to ',num2str(g.lmon.line(n).ToBus),];
% end
% 
% 
% subplot(2,1,1)
% title('Real Power Flow on Lines')
% ylabel('MW[PU]')
% xlabel('Time [seconds]')
% set(gca, 'XTick', xevents);
% legend(legNames, 'location','best')
% 
% subplot(2,1,2)
% title('Reactive Power Flow on Lines')
% ylabel('MVAR[PU]')
% xlabel('Time [seconds]')
% set(gca, 'XTick', xevents);
% legend(legNames, 'location','best')
