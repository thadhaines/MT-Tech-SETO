function livePlot(k)
%LIVEPLOT Alternate 1 - Plots fault bus v and system machine speed
% LIVEPLOT Alternate 1 - Plots fault bus v, system machine speeds,
% and lmod or pwrmod input signals if modulate detected
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

global g
if isnumeric(k)
    if (mod(k,25)==0)
        if ~isempty(g.lmod.lmod_con) || ~isempty(g.pwr.pwrmod_con)
            nPlt = 3;
        else
            nPlt = 2;
        end
        
        % plot the voltage of the faulted bus
        subplot(nPlt,1,1)
        plot(g.sys.t(1:k), abs(g.bus.bus_v(g.bus.bus_idx(1),1:k)),'r')
        busNum = int2str(g.bus.bus(g.bus.bus_idx(1),1));
        title(['Voltage Magnitude at Bus ', busNum]);
        xlabel('Time [sec]');
        ylabel('Volatge [PU]');
        
        % plot generator info
        subplot(nPlt,1,2)
        Lcolor = lines(size(g.mac.mac_spd,1));
        for pltGen = 1:size(g.mac.mac_spd,1)
            plot(g.sys.t(1:k),g.mac.mac_spd(pltGen, 1:k), 'color',Lcolor(pltGen,:))
            hold on
        end
        title('System Generator Speed');
        xlabel('Time [sec]');
        ylabel('Speed [PU]');
        
        % plot load moduation (if present)
        if ~isempty(g.lmod.lmod_con)
            subplot(nPlt,1,3)
            plot(g.sys.t(1:k), g.lmod.lmod_st(1:k))
            title('System Real Load Modulation');
            xlabel('Time [sec]');
            ylabel('MW [PU]');
        end
        
        % Plot Powermod injection if present
        if ~isempty(g.pwr.pwrmod_con)
            subplot(nPlt,1,3)
            Lcolor = lines(size(g.pwr.pwrmod_p_st,1));
            for pltData = 1:size(g.pwr.pwrmod_p_st,1)
                plot(g.sys.t(1:k),g.pwr.pwrmod_p_st(pltData, 1:k), 'color',Lcolor(pltData,:))
                hold on
            end
            title('Power Modulation');
            xlabel('Time [sec]');
            ylabel('MW [PU]');
        end
        
        drawnow
    end
    
else % end of sim full data plot
    if ~isempty(g.lmod.lmod_con) || ~isempty(g.pwr.pwrmod_con)
        nPlt = 3;
    else
        nPlt = 2;
    end
    
    % plot the voltage of the faulted bus
    subplot(nPlt,1,1)
    plot(g.sys.t, abs(g.bus.bus_v(g.bus.bus_idx(1),:)),'r')
    busNum = int2str(g.bus.bus(g.bus.bus_idx(1),1));
    title(['Voltage Magnitude at Bus ', busNum]);
    xlabel('Time [sec]');
    ylabel('Volatge [PU]');
    
    % plot generator info
    subplot(nPlt,1,2)
    Lcolor = lines(size(g.mac.mac_spd,1));
    legNames = {};
    for pltGen = 1:size(g.mac.mac_spd,1)
        plot(g.sys.t,g.mac.mac_spd(pltGen,:), 'color',Lcolor(pltGen,:))
        legNames{end+1} = ['Generator ',int2str(pltGen)];
        hold on
    end
    title('System Generator Speed');
    xlabel('Time [sec]');
    ylabel('Speed [PU]');
    legend(legNames, 'location','best')
    
    % plot load moduation (if present)
    if ~isempty(g.lmod.lmod_con)
        subplot(nPlt,1,3)
        plot(g.sys.t, g.lmod.lmod_st)
        title('System Real Load Modulation');
        xlabel('Time [sec]');
        ylabel('MW [PU]');
    end
    
    % Plot Powermod injection if present
    if ~isempty(g.pwr.pwrmod_con)
        subplot(nPlt,1,3)
        Lcolor = lines(size(g.pwr.pwrmod_p_st,1));
        for pltData = 1:size(g.pwr.pwrmod_p_st,1)
            plot(g.sys.t,g.pwr.pwrmod_p_st(pltData, :), 'color',Lcolor(pltData,:))
            hold on
        end
        title('Power Modulation');
        xlabel('Time [sec]');
        ylabel('MW [PU]');
    end
    
end
end % end function