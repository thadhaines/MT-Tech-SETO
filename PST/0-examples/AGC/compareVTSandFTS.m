% script to visually compare variable and fixed data....
% note:  For AGC

clear; close all
printFigs = 0;


caseCell = { 'AGCtestVTS.mat'};

for test=1:length(caseCell);
    
    caseN = caseCell{test};
    odeN = caseN(4:end-4);
    %%
    % load varible data
    load(caseN)
    gv = g;
    clear g;
    
    % load fixed data
    %load('FTSpstSETOdatane_hiskens.mat')
    load('AGCtestFTS.mat') % fixed - using OG s_simu_batch
    
    %% fig 1
    % bus v, bus angle
    % speed, freq
    % pss, efd
    
    figure
    %% bus voltage
    subplot(3,2,1)
    busV = 4;
    
    plot(gv.sys.t, abs(gv.bus.bus_v(busV,:)),'k','linewidth',1.25)
    hold on
    plot(g.sys.t,abs(g.bus.bus_v(busV,:)) , 'm--','linewidth',1)
    legend('VTS','FTS','location','best')
    title({ ['Bus Voltage Comparison of Bus ',int2str(busV)]})
    ylabel('Voltage [PU]')
    xlabel('Time [seconds]')
    grid on
    
    %% bus angle
    subplot(3,2,2)
    plot(gv.sys.t(1:gv.vts.dataN), unwrap(angle(gv.bus.bus_v(busV,1:gv.vts.dataN))),'k','linewidth',1.25)
    hold on
    plot(g.sys.t,unwrap(angle(g.bus.bus_v(busV,:))) , 'm--','linewidth',1)
    legend('VTS','FTS','location','best')
    title({ ['Bus Voltage Angle Comparison of Bus ',int2str(busV)]})
    ylabel('Angle [rads]')
    xlabel('Time [seconds]')
    grid on
    
    %% Machine dynamic
    subplot(3,2,3)
    mac = 1;
    
    plot(gv.sys.t, gv.mac.mac_spd(mac,:),'k','linewidth',1.25)
    hold on
    plot(g.sys.t, g.mac.mac_spd(mac,:), 'm--','linewidth',1)
    legend('VTS','FTS','location','best')
    title({ ['Machine Speed comparison of Machine ',int2str(mac)]})
    ylabel('Speed [PU]')
    xlabel('Time [seconds]')
    grid on
    
    %% system frequency average
    subplot(3,2,4)
    plot(gv.sys.t, gv.sys.aveF,'k','linewidth',1.25)
    hold on
    plot(g.sys.t, g.sys.aveF, 'm--','linewidth',1)
    legend('VTS','FTS','location','best')
    title({ ['Average System Frequency Comparison'];})
    ylabel('Frequency [PU]')
    xlabel('Time [sec]')
    grid on
    
    %% exciter dynamic
    subplot(3,2,5)
    exc = 1;
    plot(gv.sys.t, gv.exc.Efd(exc,:),'k','linewidth',1.25)
    hold on
    plot(g.sys.t, g.exc.Efd(exc,:), 'm--','linewidth',1)
    legend('VTS','FTS','location','best')
    ylabel('Voltage [PU]')
    xlabel('Time [seconds]')
    title({ ['Exciter Efd comparison of Exciter ',int2str(exc)];})
    grid on
    
    
    %% pss dynamic
    subplot(3,2,6)
    pss = 1;
    plot(gv.sys.t, gv.pss.pss_out(pss,:),'k','linewidth',1.25)
    hold on
    plot(g.sys.t, g.pss.pss_out(pss,:), 'm--','linewidth',1)
    legend('VTS','FTS','location','best')
    title({ ['PSS Output comparison of PSS ',int2str(mac)];})
    ylabel('Pss Out [PU]')
    xlabel('Time [seconds]')
    grid on
    %NE 39 bus 10 machine - 14 ms 3 Phase Fault
    
    mtit(odeN, 'xoff',0.0, 'yoff',0.025) % subplot master title - is custom function added to path
    
    %% Manipulate size and export
    if printFigs
        set(gcf, 'Position', [1 41 1366 652])
        set(gcf,'color','w'); % to remove border of figure
        export_fig([odeN,'comp'],'-pdf'); % to print fig
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
    text( 2.5/6*g.sys.t(end) , maxSln*.85,txtBlk)
    
    legend({'Variable TS', ['Fixed TS (', num2str(fts(10)), ' sec)'] }, 'location', 'northwest')
    title({odeN;'Time Step Comparison'})
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
    text(2.5/6*g.sys.t(end) , maxSln*.7,txtBlk)
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
    ylim([0,aveSln*5+0.5])
    xlim([min(gv.sys.t), max(gv.sys.t)])
    grid on
    title('Number of Network and Dynamic Solutions (Detail)')
    xlabel('Time [seconds]')
    ylabel('Number of Solutions')
    
    %% Manipulate size and export
    if printFigs
        set(gcf, 'Position', [1 41 1366 652])
        set(gcf,'color','w'); % to remove border of figure
        export_fig([odeN,'steps'],'-pdf'); % to print fig
    end
    fprintf('%s time: %s\n', odeN, gv.sys.ElapsedNonLinearTime)
    disp('')
    
    %% tg values
figure
plot(g.sys.t, g.tg.tg_sig(1,:))
hold on
plot(g.sys.t, g.tg.tg_sig(2,:))
plot(g.sys.t, g.tg.tg_sig(3,:))
plot(g.sys.t, g.tg.tg_sig(4,:))


plot(gv.sys.t, gv.tg.tg_sig(1,:),'--')
hold on
plot(gv.sys.t, gv.tg.tg_sig(2,:),'--')
plot(gv.sys.t, gv.tg.tg_sig(3,:),'--')
plot(gv.sys.t, gv.tg.tg_sig(4,:),'--')

grid on
ylabel('MW [PU]')
xlabel('Time [sec]')
title('Signal sent to governor Pref (tg\_sig)')
legend({'Gov 1 FTS', 'Gov 2 FTS','Gov 3 FTS', 'Gov 4 FTS', ...
    'Gov 1 VTS', 'Gov 2 VTS','Gov 3 VTS', 'Gov 4 VTS'},'location','best')

if printFigs
    set(gcf,'color','w'); % to remove border of figure
    export_fig([caseName,'f4'],'-pdf'); % to print fig
end

%%
figure

% % RACE
% plot(g.sys.t, g.agc.agc(1).race)
% hold on
% plot(g.sys.t, g.agc.agc(2).race)

% SACE
plot(g.sys.t, g.agc.agc(1).sace,'linewidth',2)
hold on
plot(g.sys.t, g.agc.agc(2).sace,'linewidth',2)

% % RACE
% plot(gv.sys.t, gv.agc.agc(1).race, '--')
% plot(gv.sys.t, gv.agc.agc(2).race, '--')

% SACE
plot(gv.sys.t, gv.agc.agc(1).sace, '--','linewidth',2)
plot(gv.sys.t, gv.agc.agc(2).sace, '--','linewidth',2)
legend({'SACE 1 FTS', 'SACE 2 FTS', 'SACE 1 VTS', 'SACE 2 VTS'})

%% ace sig
figure 

plot(g.sys.t, g.agc.agc(1).aceSig,'linewidth',2)
hold on
plot(g.sys.t, g.agc.agc(2).aceSig,'linewidth',2)

plot(gv.sys.t, gv.agc.agc(1).aceSig, '--','linewidth',2)
hold on
plot(gv.sys.t, gv.agc.agc(2).aceSig, '--','linewidth',2)

title('aceSig comparison')
legend({'FTS1', 'FTS2', 'VTS1', 'VTS2'})
grid on

%% ace2dist
figure 
plot(g.sys.t, g.agc.agc(1).ace2dist,'linewidth',2)
hold on
plot(g.sys.t, g.agc.agc(2).ace2dist,'linewidth',2)

plot(gv.sys.t, gv.agc.agc(1).ace2dist, '--','linewidth',2)
hold on
plot(gv.sys.t, gv.agc.agc(2).ace2dist, '--','linewidth',2)

title('ace2dist comparison')
legend({'FTS1', 'FTS2', 'VTS1', 'VTS2'})
grid on

%% ctrlgen input
figure
hold on

for n = 1:g.agc.n_agc
    for cg=1:g.agc.agc(n).n_ctrlGen
        plot(g.sys.t, g.agc.agc(n).ctrlGen(cg).input)
        plot(gv.sys.t, gv.agc.agc(n).ctrlGen(cg).input, '--')
    end
end
title('filter input')

%% ctrlgen x
figure
hold on

for n = 1:g.agc.n_agc
    for cg=1:g.agc.agc(n).n_ctrlGen
        plot(g.sys.t, g.agc.agc(n).ctrlGen(cg).x)
        plot(gv.sys.t, gv.agc.agc(n).ctrlGen(cg).x, '--')
    end
end
title('filter x')

%% ctrlgen dx
figure
hold on

for n = 1:g.agc.n_agc
    for cg=1:g.agc.agc(n).n_ctrlGen
        plot(g.sys.t, g.agc.agc(n).ctrlGen(cg).dx)
        plot(gv.sys.t, gv.agc.agc(n).ctrlGen(cg).dx, '--')
    end
end
title('filter dx')

fprintf('fixed time: %s\n', g.sys.ElapsedNonLinearTime)
end