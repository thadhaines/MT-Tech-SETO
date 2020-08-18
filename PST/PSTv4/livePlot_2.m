function livePlot(k)
%LIVEPLOT Alternate 2 - plot AGC/ACE values
% LIVEPLOT Alternate 2 - plot average system frequencies, IC error, RACE,
% SACE and DACE.
%
% Syntax: livePlot(k)
%
%   NOTES: Meant to replace livePlot function for specific case runs.
%
%   Input:
%   k - integer time (data index)
%
%   Output:
%   VOID
%
%   History:
%   Date        Time    Engineer        Description
%   08/13/20    11:14   Thad Haines     Version 1.0.0
%   08/18/20    11:21   Thad Haines     Added 'cla' command to live updates

global g
if isnumeric(k)
    if (mod(k,25)==0)
        
        % plot all area and system average frequency
        subplot(3,1,1)
        cla
        plot(g.sys.t(1:k),g.sys.aveF(1:k),'k','linewidth',2)
        if g.area.n_area > 0
            hold on
            Lcolor = lines(g.area.n_area);
            for ndx=1:g.area.n_area
                plot(g.sys.t(1:k),g.area.area(ndx).aveF(1:k),'color',Lcolor(ndx,:))
            end
        end
        title('Weighted Average Frequency');
        xlabel('Time [sec]');
        ylabel('Frequency [PU]');
        grid on
        
        % plot interchange error and RACE
        subplot(3,1,2)
        cla
        Lcolor = lines(g.agc.n_agc*2);
        for ndx = 1:g.agc.n_agc
            plot(g.sys.t(1:k),real(g.area.area(ndx).icA(1:k)) - real(g.area.area(ndx).icA(1)),'.-','linewidth',2, 'color',Lcolor(ndx,:))
            hold on
            plot(g.sys.t(1:k),g.agc.agc(ndx).race(1:k), 'color',Lcolor(ndx,:), 'color',Lcolor(ndx+1,:))
        end
        title('Area Interchange Error(thick) and AGC RACE(thin)');
        xlabel('Time [sec]');
        ylabel('MW [PU]');
        grid on
        
        % plot SACE and DACE
        subplot(3,1,3)
        cla
        Lcolor = lines(g.agc.n_agc*2);
        for ndx = 1:g.agc.n_agc
            plot(g.sys.t(1:k), g.agc.agc(ndx).sace(1:k),':','linewidth',2, 'color',Lcolor(ndx,:))
            hold on
            plot(g.sys.t(1:k), g.agc.agc(ndx).ace2dist(1:k),'--', 'color',Lcolor(ndx+1,:))
        end
        title('AGC SACE(dotted) and DACE(dashed)');
        xlabel('Time [sec]');
        ylabel('MW [PU]');
        grid on
        
        drawnow
    end
    
else % end of sim full data plot
    % same as above, but fore entire simulation time and with legends...
    
    % plot all area and system average frequency
    subplot(3,1,1)
    cla
    plot(g.sys.t,g.sys.aveF,'k','linewidth',2)
    legNames = {'System'};
    if g.area.n_area > 0
        hold on
        Lcolor = lines(g.area.n_area);
        for ndx=1:g.area.n_area
            plot(g.sys.t,g.area.area(ndx).aveF,'color',Lcolor(ndx,:))
            legNames(end+1) = {['Area ',int2str(g.area.area(ndx).number)]};
        end
    end
    title('Weighted Average Frequency');
    xlabel('Time [sec]');
    ylabel('Frequency [PU]');
    legend(legNames,'location','best')
    grid on
    
    % plot interchange error and RACE
    subplot(3,1,2)
    cla
    Lcolor = lines(g.agc.n_agc*2);
    legNames = {};
    for ndx = 1:g.agc.n_agc
        plot(g.sys.t,real(g.area.area(ndx).icA) - real(g.area.area(ndx).icA(1)),'.-','linewidth',2, 'color',Lcolor(ndx,:))
        hold on
        plot(g.sys.t, g.agc.agc(ndx).race, 'color',Lcolor(ndx,:), 'color',Lcolor(ndx+1,:))
        legNames(end+1) = {['Area ',int2str(g.area.area(ndx).number),' IC Error']};
        legNames(end+1) = {['Area ',int2str(g.area.area(ndx).number),' RACE']};
    end
    title('Area Interchange Error(thick) and AGC RACE(thin)');
    xlabel('Time [sec]');
    ylabel('MW [PU]');
    legend(legNames,'location','best')
    grid on
    
    % plot SACE and DACE
    subplot(3,1,3)
    cla
    Lcolor = lines(g.agc.n_agc*2);
    legNames={};
    for ndx = 1:g.agc.n_agc
        plot(g.sys.t, g.agc.agc(ndx).sace,':','linewidth',2, 'color',Lcolor(ndx,:))
        hold on
        plot(g.sys.t, g.agc.agc(ndx).ace2dist,'--', 'color',Lcolor(ndx+1,:))
        legNames(end+1) = {['Area ',int2str(g.area.area(ndx).number),' SACE']};
        legNames(end+1) = {['Area ',int2str(g.area.area(ndx).number),' DACE']};
    end
    title('AGC SACE(dotted) and DACE(dashed)');
    xlabel('Time [sec]');
    ylabel('MW [PU]');
    legend(legNames,'location','best')
    grid on
    
    drawnow
    
end
end % end function