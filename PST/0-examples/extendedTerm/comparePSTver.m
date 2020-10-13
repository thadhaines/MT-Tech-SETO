% script to visually compare variable and fixed data....
% note:  For extended term sim

clear; close all
printFigs = 0;
detailPlots = 1;
detailX = [54,64]; % detail 1
detailX = [119,129]; % detail 2

% load varible data
% load('AGCtestVTS.mat')
load('pmRampVTS.mat') % original tolerances Of 'RelTol',1e-4,'AbsTol',1e-7 -61 seconds
%load('pmRampFTS.mat') % increased 'RelTol',1e-5,'AbsTol',1e-7, - 106 sec
% load('pmRampVTS2.mat') % altered time blocks
gv = g;

% load seto data
load('pmRampSETO.mat') % fixed
clearvars -except gv g printFigs detailPlots detailX

% load 311 data
load pmRamp3FTS.mat bus_v mac_spd Efd pss_out t pmech pm_sig pelect ets

%% fig 1
% bus v, bus angle
% speed, freq
% pss, efd

figure
% bus voltage
subplot(3,2,1)
busV = 4;

plot(g.sys.t, abs(g.bus.bus_v(busV,:)),'k','linewidth',1.25)
hold on
plot(gv.sys.t,abs(gv.bus.bus_v(busV,:)) , 'm--','linewidth',1)
plot(t,abs(bus_v(busV,:)) , 'c-.','linewidth',1)

legend('SETO','VTS','3.1','location','best')
title({ ['Bus Voltage Comparison of Bus ',int2str(busV)]})
ylabel('Voltage [PU]')
xlabel('Time [seconds]')
grid on

if detailPlots
    xlim(detailX)
end
% bus angle
% locate slack bus angle
sbAng = g.bus.bus( g.bus.bus(:,10) == 1);

subplot(3,2,2)
plot(g.sys.t, rad2deg( unwrap(angle(g.bus.bus_v(busV,:))) - unwrap(angle(g.bus.bus_v(sbAng,:)))) ,'k','linewidth',1.25)
hold on
plot(gv.sys.t,rad2deg( unwrap(angle(gv.bus.bus_v(busV,:))) - unwrap(angle(gv.bus.bus_v(sbAng,:)))) , 'm--','linewidth',1)
plot(t,rad2deg( unwrap(angle(bus_v(busV,:))) - unwrap(angle(bus_v(sbAng,:)))) , 'c-.','linewidth',1)
legend('SETO','VTS','3.1','location','best')
title({ ['Bus Voltage Angle Comparison of Bus ',int2str(busV)]})
ylabel('Angle [degree]')
xlabel('Time [seconds]')
grid on
if detailPlots
    xlim(detailX)
end

%% Machine dynamic
subplot(3,2,3)
mac = 1;

plot(g.sys.t, g.mac.mac_spd(mac,:),'k','linewidth',1.25)
hold on
plot(gv.sys.t, gv.mac.mac_spd(mac,:), 'm--','linewidth',1)
plot(t, mac_spd(mac,:), 'c-.','linewidth',1)
legend('SETO','VTS','3.1','location','best')
title({ ['Machine Speed comparison of Machine ',int2str(mac)]})
ylabel('Speed [PU]')
xlabel('Time [seconds]')
grid on
if detailPlots
    xlim(detailX)
end

%% system frequency average
subplot(3,2,4)
plot(g.sys.t, g.sys.aveF,'k','linewidth',1.25)
hold on
plot(gv.sys.t, gv.sys.aveF, 'm--','linewidth',1)

legend('SETO','VTS','location','best')
title({ ['Average System Frequency Comparison'];})
ylabel('Frequency [PU]')
xlabel('Time [sec]')
grid on
if detailPlots
    xlim(detailX)
end

%% exciter dynamic
subplot(3,2,5)
exc = 1;
plot(g.sys.t, g.exc.Efd(exc,:),'k','linewidth',1.25)
hold on
plot(gv.sys.t, gv.exc.Efd(exc,:), 'm--','linewidth',1)
plot(t, Efd(exc,:), 'c-.','linewidth',1)

legend('SETO','VTS','3.1','location','best')
ylabel('Voltage [PU]')
xlabel('Time [seconds]')
title({ ['Exciter Efd comparison of Exciter ',int2str(exc)];})
grid on
if detailPlots
    xlim(detailX)
end

%% pss dynamic
subplot(3,2,6)
pss = 1;
plot(g.sys.t, g.pss.pss_out(pss,:),'k','linewidth',1.25)
hold on
plot(gv.sys.t, gv.pss.pss_out(pss,:), 'm--','linewidth',1)
plot(t, pss_out(pss,:), 'c-.', 'linewidth', 1)

legend('SETO','VTS', '3.1', 'location','best')
title({ ['PSS Output comparison of PSS ',int2str(mac)];})
ylabel('Pss Out [PU]')
xlabel('Time [seconds]')
grid on
if detailPlots
    xlim(detailX)
end
%mtit(odeN, 'xoff',0.0, 'yoff',0.025) % subplot master title - is custom function added to path

%% Manipulate size and export
if printFigs
    set(gcf, 'Position', [1 41 1366 652])
    set(gcf,'color','w'); % to remove border of figure
    export_fig(['PWRcomp'],'-pdf'); % to print fig
end

%% time step size calculation
vts = zeros(gv.vts.dataN-1,1);
for n=2:gv.vts.dataN
    vts(n-1)= gv.sys.t(n)-gv.sys.t(n-1);
end

fts = zeros(size(g.sys.t,2)-1,1);
for n=2:size(g.sys.t,2)
    fts(n-1) = g.sys.t(n) - g.sys.t(n-1);
end

% time step plottings
figure
subplot(3,1,1)
stairs(gv.sys.t(1:end-1),vts,'k','linewidth',1.25)
hold on
grid on
stairs(g.sys.t(1:end-1), fts, 'm--','linewidth',1)


aveSln = mean(vts);
maxSln = max(vts);

minSln = aveSln;
for i=1:length(vts)
    if vts(i) ~= 0
        if vts(i) < minSln
            minSln = vts(i);
        end
    end
end


txtBlk = {['Average Time Step:  ', num2str(aveSln,3)]; ...
    ['Maximum Time Step: ', num2str(maxSln,3)];  ...
    ['Minimum Time Step: ', num2str(minSln,3)]   };
grid on
text( 150, 0.15, txtBlk)
xlim([0, 240])

legend({'Variable TS', ['Fixed TS (', num2str(fts(10)), ' sec)'] }, 'location', 'northwest')
title({'Time Step Comparison'})
ylabel('Time Step Size [seconds]')
xlabel('Time [seconds]')
ylim([0, .3])
xlim([0, 240])

% plot number of solutions
subplot(3,1,2)
stairs(gv.sys.t, gv.vts.slns,'k','linewidth',1.25)
hold on
bar(gv.sys.t, gv.vts.slns,'k')

aveSln = round(mean(gv.vts.slns));
maxSln = max(gv.vts.slns);
txtBlk = {['Average Solutions per step: ', int2str(aveSln)]; ...
    ['Maximum Solutions per step: ', int2str(maxSln)]; ...
    ['Total Solutions: ', int2str(gv.vts.tot_iter)]; ...
    ['Total Steps: ', int2str(gv.vts.dataN)]         };
grid on

text( 150, 60, txtBlk)
xlim([0, 240])
title('Number of Network and Dynamic Solutions')
%xlabel('Time [seconds]')
%ylim([aveSln*10, maxSln*1.2])
ylabel('Number of Solutions')

xlim([0, 240])
% plot number of solutions detail
subplot(3,1,3)
stairs(gv.sys.t, gv.vts.slns,'k','linewidth',1.25)
hold on
bar(gv.sys.t, gv.vts.slns,'k')
ylim([0,aveSln*2+0.5])
xlim([min(gv.sys.t), max(gv.sys.t)])
grid on
title('Number of Network and Dynamic Solutions (Detail)')
xlabel('Time [seconds]')
ylabel('Number of Solutions')
xlim([0, 240])

%% Manipulate size and export
if printFigs
    set(gcf, 'Position', [1 41 1366 652])
    set(gcf,'color','w'); % to remove border of figure
    export_fig(['PWRsteps'],'-pdf'); % to print fig
end
%%
fprintf('VTS time: %s\n', gv.sys.ElapsedNonLinearTime)
fprintf('SETO time: %s\n', g.sys.ElapsedNonLinearTime)
fprintf('3.1 time: %s\n', ets)

%% mechanical power
figure
legNames = {};
for n=1:size(g.mac.pmech,1)
plot(g.sys.t, g.mac.pmech(n,:))
hold on
plot(gv.sys.t, gv.mac.pmech(n,:) ,'--')
plot(t, pmech(n,:) ,'-.')
legNames = [legNames, ['SETO gen ',int2str(n)], ['VTS gen ',int2str(n)], ['3.1.1 gen ',int2str(n)]];
end

grid on
xlabel('Time [seconds]')
ylabel('MW [PU]')
title('Mechanical Power')
legend(legNames,'location','best')

%% pm sig
figure
legNames={};
for n=1:size(g.mac.pm_sig,1)
plot(g.sys.t, g.mac.pm_sig(n,:))
hold on
plot(gv.sys.t, gv.mac.pm_sig(n,:) ,'--')
plot(t, pm_sig(n,:) ,'-.')
legNames = [legNames, ['SETO gen ',int2str(n)], ['VTS gen ',int2str(n)], ['3.1.1 gen ',int2str(n)]];
end
grid on
xlabel('Time [seconds]')
ylabel('MW [PU]')
title('Mechanical Power Modulation Signal')
legend(legNames,'location','best')

%% p elect
figure
legNames={};
for n=1:size(g.mac.pelect,1)
plot(g.sys.t, g.mac.pelect(n,:))
hold on
plot(gv.sys.t, gv.mac.pelect(n,:) ,'--')
plot(t, pelect(n,:) ,'-.')
legNames = [legNames, ['SETO gen ',int2str(n)], ['VTS gen ',int2str(n)], ['3.1.1 gen ',int2str(n)]];
end
grid on
xlabel('Time [seconds]')
ylabel('MW [PU]')
title('Electric Power')
legend(legNames,'location','best')
% %%
% figure
% plot(t, pmech)
% 
% figure
% plot(t, pm_sig)
% 
% figure
% plot(t, pelect)
% grid on

ylabel('MW [PU]')
xlabel('Time [sec]')

if printFigs
    set(gcf,'color','w'); % to remove border of figure
    export_fig(['pwrSig'],'-pdf'); % to print fig
end


