% compare 10 miute AGC miniWECC recovery
clear; close all

%load('VTSagcPstep') % relative tolerance to 1e-4
%load('VTSrel-5') % decreased relative tolerance to 1e-5
load('test') % decreased relative tolerance to 1e-5, decreased number of time blocks...
gv = g;
clearvars -except gv

load('FTSagcPstep')

clearvars -except gv g
%%
fprintf('Elapsed Time (Huens):\t %s\n', g.sys.ElapsedNonLinearTime)
fprintf('Elapsed Time (ode23):\t %s\n', gv.sys.ElapsedNonLinearTime)
fprintf('VTS speed up: %4.4f\n', ...
    str2double(g.sys.ElapsedNonLinearTime)/str2double(gv.sys.ElapsedNonLinearTime))
fprintf('VTS logged data reduction: %4.4f\n', ...
    1172615366/32419352) % collected from whos('varName')


%% copied time step plot
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
    text(100, maxSln*.7,txtBlk)
    
    legend({'Variable TS', ['Fixed TS (', num2str(fts(10)), ' sec)'] }, 'location', 'north')
    title({'Time Step Comparison'})
    ylabel('Time Step Size [seconds]')
    xlabel('Time [seconds]')
    
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
    text(100, maxSln*.7,txtBlk)
    title('Number of Network and Dynamic Solutions')
    %xlabel('Time [seconds]')
    %ylim([aveSln*10, maxSln*1.2])
    ylabel('Number of Solutions')
    
    xlim([min(gv.sys.t), max(gv.sys.t)])
    
    % plot number of solutions detail
    subplot(3,1,3)
    stairs(gv.sys.t, gv.vts.slns,'k','linewidth',1.25)
    hold on
    bar(gv.sys.t, gv.vts.slns,'k')
    ylim([0,aveSln+0.5])
    xlim([min(gv.sys.t), max(gv.sys.t)])
    grid on
    title('Number of Network and Dynamic Solutions (Detail)')
    xlabel('Time [seconds]')
    ylabel('Number of Solutions')
    
%% multiple comparison figure
plotLims = [0,600];
figure
% bus voltage
subplot(3,2,1)
busV = 2;

plot(g.sys.t, abs(g.bus.bus_v(busV,:)),'k','linewidth',1.25)
hold on
plot(gv.sys.t,abs(gv.bus.bus_v(busV,:)) , 'm--','linewidth',1)
legend('FTS','VTS','location','best')
title({ ['Bus Voltage Comparison of Bus ',int2str(busV)]})
ylabel('Voltage [PU]')
xlabel('Time [seconds]')
xlim(plotLims)
grid on

% bus angle
% locate slack bus angle
sbAng = g.bus.bus( g.bus.bus(:,10) == 1);

subplot(3,2,2)
plot(g.sys.t, rad2deg( unwrap(angle(g.bus.bus_v(busV,:))) - unwrap(angle(g.bus.bus_v(sbAng,:)))) ,'k','linewidth',1.25)
hold on
plot(gv.sys.t,rad2deg( unwrap(angle(gv.bus.bus_v(busV,:))) - unwrap(angle(gv.bus.bus_v(sbAng,:)))) , 'm--','linewidth',1)
legend('FTS','VTS','location','best')
title({ ['Bus Voltage Angle Comparison of Bus ',int2str(busV)]})
ylabel('Angle [degree]')
xlabel('Time [seconds]')
xlim(plotLims)
grid on

% Machine dynamic
subplot(3,2,3)
mac = 1;

plot(g.sys.t, g.mac.mac_spd(mac,:),'k','linewidth',1.25)
hold on
plot(gv.sys.t, gv.mac.mac_spd(mac,:), 'm--','linewidth',1)
legend('FTS','VTS','location','best')
title({ ['Machine Speed comparison of Machine ',int2str(mac)]})
ylabel('Speed [PU]')
xlabel('Time [seconds]')
xlim(plotLims)
grid on

% system frequency average
subplot(3,2,4)
plot(g.sys.t, g.sys.aveF,'k','linewidth',1.25)
hold on
plot(gv.sys.t, gv.sys.aveF, 'm--','linewidth',1)
legend('FTS','VTS','location','best')
title({ ['Average System Frequency Comparison'];})
ylabel('Frequency [PU]')
xlabel('Time [sec]')
xlim(plotLims)
grid on

% exciter dynamic
subplot(3,2,5)
exc = 1;
plot(g.sys.t, g.exc.Efd(exc,:),'k','linewidth',1.25)
hold on
plot(gv.sys.t, gv.exc.Efd(exc,:), 'm--','linewidth',1)
legend('FTS','VTS','location','best')
ylabel('Voltage [PU]')
xlabel('Time [seconds]')
title({ ['Exciter Efd comparison of Exciter ',int2str(exc)];})
xlim(plotLims)
grid on


% pss dynamic
subplot(3,2,6)
pss = 1;
plot(g.sys.t, g.pss.pss_out(pss,:),'k','linewidth',1.25)
hold on
plot(gv.sys.t, gv.pss.pss_out(pss,:), 'm--','linewidth',1)
legend('FTS','VTS','location','best')
title({ ['PSS Output comparison of PSS ',int2str(pss)];})
ylabel('Pss Out [PU]')
xlabel('Time [seconds]')
xlim(plotLims)
grid on

%%
figure
subplot(3,1,1)

for n=1:3
    plot(g.sys.t, g.agc.agc(n).race, 'linewidth',2)
    hold on
    plot(gv.sys.t, gv.agc.agc(n).race, '--', 'linewidth',1.5)
end
title('RACE Comparisons')
grid on
ylabel('MW [PU]')
xlabel('Time [sec]')

subplot(3,1,2)

for n=1:3
    plot(g.sys.t, real(g.area.area(n).icA-g.area.area(n).icS), 'linewidth',2)
    hold on
    plot(gv.sys.t, real(gv.area.area(n).icA-gv.area.area(n).icS), '--', 'linewidth',1.5)
end
grid on
title('Area Interchange Comparisons')
ylabel('MW [PU]')
xlabel('Time [sec]')

subplot(3,1,3)
for n = 1:g.agc.n_agc
    for cg=1:g.agc.agc(n).n_ctrlGen
        plot(g.sys.t, g.agc.agc(n).ctrlGen(cg).input)
        hold on
        plot(gv.sys.t, gv.agc.agc(n).ctrlGen(cg).input, '--')
    end
end
title('Controlled Machine Comparisons')
grid on
ylabel('MW [PU]')
xlabel('Time [sec]')