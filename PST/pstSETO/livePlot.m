%% fancier live plotting ~ Obviously causes sim to run slower

if (mod(k,25)==0)
    if ~isempty(g.lmod.lmod_con) || ~isempty(pwrmod_con)
        nPlt = 3;
    else
        nPlt = 2;
    end
    
    % format bus voltage for plot
    v_p(1:k)=abs(bus_v(bus_idx(1),1:k));
    % plot the voltage of the faulted bus
    subplot(nPlt,1,1)
    plot(t(1:k),v_p(1:k),'r')
    title('Voltage Magnitude at Fault Bus');
    xlabel('Time [sec]');
    ylabel('Volatge [PU]');
    
    % plot generator info
    subplot(nPlt,1,2)
    Lcolor = lines(size(g.mac.mac_spd,1));
    for pltGen = 1:size(g.mac.mac_spd,1)
        plot(t(1:k),g.mac.mac_spd(pltGen, 1:k), 'color',Lcolor(pltGen,:))
        hold on
    end
    title('System Generator Speed');
    xlabel('Time [sec]');
    ylabel('Speed [PU]');
    
    % plot load moduation (if present)
    if ~isempty(g.lmod.lmod_con)
        subplot(nPlt,1,3)
        plot(t(1:k),g.lmod.lmod_st(1:k))
        title('System Real Load Modulation');
        xlabel('Time [sec]');
        ylabel('MW [PU]');
    end
    
    % Plot Powermod injection if present
    if ~isempty(pwrmod_con)
        subplot(nPlt,1,3)
        Lcolor = lines(size(pwrmod_p_st,1));
        for pltData = 1:size(pwrmod_p_st,1)
            plot(t(1:k),pwrmod_p_st(pltData, 1:k), 'color',Lcolor(pltData,:))
            hold on
        end
        title('Power Modulation');
        xlabel('Time [sec]');
        ylabel('MW [PU]');
    end
    
    drawnow
end