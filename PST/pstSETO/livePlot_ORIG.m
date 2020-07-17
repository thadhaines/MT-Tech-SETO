function livePlot(k)
%% live plotting ~ Obviously causes sim to run slower
global g
if isnumeric(k)
    if (mod(k,25)==0)
        % plot the voltage of the faulted bus
        plot(g.sys.t(1:k), abs(g.sys.bus_v(g.bus.bus_idx(1),1:k)),'r')
        busNum = int2str(g.bus.bus(g.bus.bus_idx(1),1));
        title(['Voltage Magnitude at Bus ', busNum]);
        xlabel('Time [sec]');
        ylabel('Volatge [PU]');
        
        drawnow
    end
else
    % plot the voltage of the faulted bus for all time - executed at end of sim
    plot(g.sys.t, abs(g.sys.bus_v(g.bus.bus_idx(1),:)),'r')
    busNum = int2str(g.bus.bus(g.bus.bus_idx(1),1));
    title(['Voltage Magnitude at Bus ', busNum]);
    xlabel('Time [sec]');
    ylabel('Volatge [PU]');
    
    drawnow
    
end