% Run 3mac_untripTest
% Only functional in PSTv4

clear all; close all; clc

%% Add correct PST verstion path to MATLAB in a generic way
folderDepth = 2; % depth of current directory from main PST directory
pstVer =   'PSTv4'; % 'pstV3p1'; % 'pstV2P3'; % 'pstSETO';%
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

copyfile('d_3mac_untripTest.m',[PSTpath 'DataFile.m']); % system data file
copyfile([PSTpath 'livePlot_1.m'],[PSTpath 'livePlot.m']); % set live plot

% % Move correct load modulation
% if strcmp(pstVer, 'pstSETO') || strcmp(pstVer, 'PSTv4')
%     copyfile('ml_sig_smallStepG.m', [PSTpath 'ml_sig.m']);
% else
%     error('not handled yet') %use of individual globals
% end

% Move correct trip logic files
if strcmp(pstVer, 'pstSETO') || strcmp(pstVer, 'PSTv4')
    copyfile('mac_trip_logic_Gen_3_G.m', [PSTpath 'mac_trip_logic.m']);
    copyfile('mtg_sig_PrefRamp.m', [PSTpath 'mtg_sig.m']);
    %copyfile('mpm_sig_PmRampG.m', [PSTpath 'mpm_sig.m']);
else
    error('not handled yet')
end

livePlotFlag = 1; % not ideal for extended term sim - may cause crashes/freezing

if strcmp(pstVer , 'PSTv4')
    s_simu
else
    s_simu_Batch %Run PST
end

%% Clean up modulation file alterations.
copyfile([PSTpath 'mac_trip_logic_ORIG.m'], [PSTpath 'mac_trip_logic.m']);
copyfile([PSTpath 'livePlot_ORIG.m'],[PSTpath 'livePlot.m']);
copyfile([PSTpath 'ml_sig_ORIG.m'],[PSTpath 'ml_sig.m']);
copyfile([PSTpath 'mac_sub_ORIG.m'],[PSTpath 'mac_sub.m']);
copyfile([PSTpath 'mpm_sig_ORIG.m'], [PSTpath 'mpm_sig.m']);
copyfile([PSTpath 'mtg_sig_ORIG.m'], [PSTpath 'mtg_sig.m']);

%% Save cleaned output data
save('3mac_untripTest_noEXC.mat'); %Save simulation outputs

%% temp file clean up
delete('PSTpath.mat')

%% ===============================================================

%% Plots
xevents = g.sys.sw_con(:,1)';
plotCell = { ...
    %f1, f2
    'mac','pmech';
    'tg','tg_sig';
    'mac','mac_spd';
    'mac','pelect';
    'mac','qelect';
    %'mac','cur_re';
    %'mac','cur_im';
    %'mac','mac_ang';
    %'mac','psi_re';
    %'mac','psi_im';
    };

% nS = find(g.sys.t > 24);
% nS = nS(1);
% nE = find(g.sys.t > 28);
% nE = nE(1);

axLim = [24,28];

lnClr=[0,0,0; 0.66, 0.66, 0.66; 1, 0, 1]; % for custom colors

for n = 1:size(plotCell,1)
    f1 = plotCell{n,1};
    f2 = plotCell{n,2};
    figure
    
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

%
% %% tripped gen xfmr bus voltages
% figure
% plot(g.sys.t, abs(g.bus.bus_v(6,:)))
% hold on
% plot(g.sys.t, abs(g.bus.bus_v(7,:)), '--')
% title('Generator and XFMR Bus Voltages')
% ylabel('Bus Voltage [PU]')
% xlabel('Time [seconds]')
% legend({'High Side XFRMR','Gen Side XFRMR'},'location','best')
% grid on

%% System Bus voltages
figure
legNames = {};
for n=1:size(g.bus.bus_v,1)-1
    plot(g.sys.t, abs(g.bus.bus_v(n,:)))
    legNames{end+1} = ['Bus ',num2str(g.bus.bus(n,1))];
    hold on
end

title('All Sustem Bus Voltages')
ylabel('Bus Voltage [PU]')
xlabel('Time [seconds]')
set(gca, 'XTick', xevents);
legend(legNames, 'location','best')
grid on

%% Line currents feeding load
figure
legNames = {};
for n=1:g.lmon.n_lmon
    subplot(2,1,1)
    plot(g.sys.t, real(g.lmon.line(n).sFrom),'color',lnClr(n,:),'linewidth',1.25)
    hold on
    grid on
    
    subplot(2,1,2)
    plot(g.sys.t, imag(g.lmon.line(n).sFrom),'color',lnClr(n,:),'linewidth',1.25)
    hold on
    grid on
    
    legNames{end+1} = ['Line ',num2str(g.lmon.line(n).FromBus), ' to ',num2str(g.lmon.line(n).ToBus),];
end


subplot(2,1,1)
title('Real Power Flow on Lines')
ylabel('MW[PU]')
xlabel('Time [seconds]')
set(gca, 'XTick', xevents);
legend(legNames, 'location','best')

subplot(2,1,2)
title('Reactive Power Flow on Lines')
ylabel('MVAR[PU]')
xlabel('Time [seconds]')
set(gca, 'XTick', xevents);
legend(legNames, 'location','best')
