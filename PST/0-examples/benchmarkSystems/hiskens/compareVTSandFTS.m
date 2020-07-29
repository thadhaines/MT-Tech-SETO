% script to visually compare variable and fixed data....
% note:  a better naming scheme could prevent numerous loads of fixed data...

clear; close all
printFigs = 0;

caseCell = { 'VTSode113.mat', 'VTSode15s.mat', 'VTSode23.mat', ...
    'VTSode23t.mat', 'VTSode23tb.mat'};

for test=1:length(caseCell);
    
    caseN = caseCell{test};
    odeN = caseN(4:end-4);
    %%
    % load varible data
    load(caseN)
    gv = g;
    clear g;
    
    % load fixed data
    load('FTSpstSETOdatane_hiskens.mat')
    
    %% fig 1
    % bus v, bus angle
    % speed, freq
    % pss, efd
    
    figure
    %% bus voltage
    subplot(3,2,1)
    busV = 16;
    
    plot(gv.sys.t(1:gv.vts.dataN), abs(gv.bus.bus_v(busV,1:gv.vts.dataN)),'k','linewidth',1.25)
    hold on
    plot(g.sys.t,abs(g.bus.bus_v(busV,:)) , 'm--','linewidth',1)
    legend('VTS','FTS','location','best')
    title({ ['Bus Voltage Comparison of Bus ',int2str(busV)]})
    ylabel('Voltage [PU]')
    xlabel('Time [seconds]')
    grid on
    
    %% bus angle
    subplot(3,2,2)
    plot(gv.sys.t(1:gv.vts.dataN), angle(gv.bus.bus_v(busV,1:gv.vts.dataN)),'k','linewidth',1.25)
    hold on
    plot(g.sys.t,angle(g.bus.bus_v(busV,:)) , 'm--','linewidth',1)
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
    mtit(odeN, 'xoff',0.0, 'yoff',0.025)
    
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
    subplot(2,1,1)
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
    text(3, maxSln*.85,txtBlk)
    
    legend({'Variable TS', ['Fixed TS (', num2str(fts(10)), ' sec)'] }, 'location', 'north')
    title({odeN;'Time Step Comparison'})
    ylabel('Time Step Size [seconds]')
    xlabel('Time [seconds]')
    
    % plot number of solutions
    subplot(2,1,2)
    stairs(gv.sys.t, gv.vts.slns,'k','linewidth',1.25)
    
    aveSln = round(mean(gv.vts.slns));
    maxSln = max(gv.vts.slns);
    txtBlk = {['Average Solutions per step: ', int2str(aveSln)]; ...
        ['Maximum Solutions per step: ', int2str(maxSln)]; ...
        ['Total Solutions: ', int2str(gv.vts.tot_iter)]; ...
        ['Total Steps: ', int2str(gv.vts.dataN)]         };
    grid on
    text(7, maxSln*.7,txtBlk)
    title('Number of Network and Dynamic Solutions')
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
end
fprintf('fixed time: %s\n', g.sys.ElapsedNonLinearTime)