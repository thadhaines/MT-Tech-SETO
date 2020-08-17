function standAlonePlot(scriptRunFlag)
%STANDALONEPLOT performs interative ploting ala OG s_simu
% STANDALONEPLOT performs interative ploting ala OG s_simu
%
% Syntax: standAlonePlot(scriptRunFlag)
%
%   NOTES: will act if scriptRunFlag == 1
%
%   Input:
%   scriptRunFlag - sets original flag to enter while loop
%
%   Output:
%   VOID
%
%   History:
%   Date        Time    Engineer        Description
%   07/30/20    10:54   Thad Haines     Version 1

%% Remaining 'loose' globals

% DeltaP/omega filter variables - 21
global  dpw_con dpw_out dpw_pot dpw_pss_idx dpw_mb_idx dpw_idx n_dpw dpw_Td_idx dpw_Tz_idx
global  sdpw1 sdpw2 sdpw3 sdpw4 sdpw5 sdpw6
global  dsdpw1 dsdpw2 dsdpw3 dsdpw4 dsdpw5 dsdpw6

% pss design - 3 - Not used in Simulation? - thad 07/18/20
global ibus_con  netg_con  stab_con

%%
global g

if scriptRunFlag == 1
    flag = 0;
    disp('***')
    disp('*** Welcom to Ineractive Plotting!')
else
    flag = 1;
end
%%
while(flag == 0)
    disp('***')
    disp('Enter number:')
    disp('     1 to plot all machine angles in 3D')
    disp('     2 to plot all machine speed deviation in 3D')
    disp('     3 to plot all machine turbine powers')
    disp('     4 to plot all machine electrical powers')
    disp('     5 to plot all exciter field voltages')
    disp('     6 to plot all bus voltage magnitude in 3D')
    disp('     7 to plot the line power flows')
    disp('     0 to quit')
    sel = input('Enter Selection >> ');
    if isempty(sel)
        sel = 0;
    end
    if sel == 1
        figure
        mesh(g.sys.t,1:1:g.mac.n_mac,g.mac.mac_ang*180/pi)
        title('Machine Angles')
        xlabel('Time [seconds]')
        ylabel('Internal Generator Number')
        zlabel('Angle [degrees]')
    elseif sel == 2
        figure
        lt = length(g.sys.t);
        mesh(g.sys.t, 1:1:g.mac.n_mac, g.mac.mac_spd- ones(g.mac.n_mac,lt) )
        title('Machine Speed Deviations')
        xlabel('Time [seconds]')
        ylabel('Internal Generator Number')
        zlabel('Speed Deviation [PU]')
    elseif sel == 3
        figure
        plot(g.sys.t,g.mac.pmech)
        grid on
        title('Turbine Power')
        xlabel('Time [seconds]')
        ylabel('MW [PU]')
    elseif sel == 4
        figure
        plot(g.sys.t,g.mac.pelect)
        grid on
        title('Generator Electric Power')
        xlabel('Time [seconds]')
        ylabel('MW [PU]')
    elseif sel == 5
        figure
        plot(g.sys.t,g.exc.Efd)
        grid on
        title('Exciter Field Voltages')
        xlabel('Time [seconds]')
        ylabel('Voltage [PU]')
    elseif sel == 6
        figure
        nbus= size(g.bus.bus_v,1);
        mesh(g.sys.t,(1:1:nbus),abs(g.bus.bus_v))
        xlabel('Time [seconds]')
        ylabel('Bus Number')
        zlabel('Voltage [PU]')
        clear nbus
    elseif sel == 7
        figure
        nline = length(g.line.line(:,1));
        if nline<25
            V1 = g.bus.bus_v(g.bus.bus_int(g.line.line(:,1)),:);
            V2 = g.bus.bus_v(g.bus.bus_int(g.line.line(:,2)),:);
            R = g.line.line(:,3);
            X = g.line.line(:,4);
            B = g.line.line(:,5);
            tap = g.line.line(:,6);
            phi = g.line.line(:,7);
        else
            % ask for lines to be plotted
            disp('Enter a single line, or rangle of lines:')
            disp('(for example: 2 or 4:5 or [4, 7, 10] )')
            line_range = input('>> ');
            if isempty(line_range)
                line_range = 1:round(size(g.line.line,1)/8);
            end
            V1 = g.bus.bus_v(g.bus.bus_int(g.line.line(line_range,1)),:);
            V2 = g.bus.bus_v(g.bus.bus_int(g.line.line(line_range,2)),:);
            R = g.line.line(line_range,3);
            X = g.line.line(line_range,4);
            B = g.line.line(line_range,5);
            tap = g.line.line(line_range,6);
            phi = g.line.line(line_range,7);
        end
        
        [S1,S2] = line_pq(V1,V2,R,X,B,tap,phi);
        plot(g.sys.t,real(S1));
        grid on
        title('Line Real Power Flow')
        xlabel('Time [seconds]')
        ylabel('MW [PU]')
        
        clear V1 V2 R X B tap phi S1 S2
        
    elseif sel == 0
        flag = 1;
    else
        error('invalid selection...')
    end
end
disp('***')
end% end function