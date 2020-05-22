%% Plotting of PST outputs
% assumes data of interst has been loaded - will fail to plot data that doesn't exist
% intended as an example / WIP plot function sandbox

%% Synchronous generator info
figure
plot(t,pmech)
title([{'pmech'}; {'generator mechanical power output'}])
xlabel('Time [sec]')
ylabel('MW [PU]')

figure
plot(t,pelect)
title([{'pelect'}; {'generator active power output'}])
xlabel('Time [sec]')
ylabel('MW [PU]')

figure
plot(t,qelect)
title([{'qelect'}; {'generator reactive power output'}])
xlabel('Time [sec]')
ylabel('MVAR [PU]')

figure
plot(t,mac_spd)
title([{'mac\_spd'}; {'machine speed'}])
xlabel('Time [sec]')
ylabel('Speed [PU]')

figure
plot(t,mac_ang)
title([{'mac\_ang'}; {'machine angle'}])
xlabel('Time [sec]')
ylabel('Angle [rads]')

%% Exciter info
figure
plot(t,Efd)
title([{'Efd'}; {'excitation output voltage (field voltage)'}])
xlabel('Time [sec]')
ylabel('Voltage [PU]')

%% bus info
% extra bus = fault bus?
figure
plot(t,abs(bus_v))
title([{'abs(bus\_v)'}; {'bus voltage magnitude'}])
xlabel('Time [sec]')
ylabel('Voltage [PU]')

figure
plot(t,angle(bus_v))
title([{'angle(bus\_v)'}; {'bus voltage angle'}])
xlabel('Time [sec]')
ylabel('Angle [rads]')

%% Branch currents
% magnitues should be the same to and from
figure
plot(t,abs(ilf))
title([{'abs(ilf)'}; {'line flow from bus current magnitude'}])
xlabel('Time [sec]')
ylabel('Current [PU]')

figure
plot(t,abs(ilt))
title([{'abs(ilt)'}; {'line flow to bus current magnitude'}])
xlabel('Time [sec]')
ylabel('Current [PU]')

% angles should be 180 degree out of phase...
figure
plot(t,rad2deg(angle(ilf(4,:))))
title([{'angle(ilf(4,:)))'}; {'line flow from bus current angle'}])
xlabel('Time [sec]')
ylabel('Angle [rads]')

figure
plot(t,rad2deg(angle(ilt(4,:))))
title([{'angles(ilt(4,:)))'}; {'line flow to bus current angle'}])
xlabel('Time [sec]')
ylabel('Angle [rads]')

figure
plot(t,rad2deg(angle(ilt(4,:)))-rad2deg(angle(ilf(4,:))))
title([{'rad2deg(angle(ilt(4,:)))-rad2deg(angle(ilf(4,:)))'}; {'line flow current angle dif'}])
xlabel('Time [sec]')
ylabel('Angle [rads]')

%% P and Q compuation investigation
figure
plot(t, real(sInjT))
title([{'real(sInjT)'}; {'active power line flow'}])
xlabel('Time [sec]')
ylabel('MW [PU]?')

figure
plot(t, imag(sInjT))
title([{'imag(sInjT)'}; {'reactive power line flow'}])
xlabel('Time [sec]')
ylabel('MVAR [PU]?')

%% Modulated real power loads

if exist('lmod_con','var')
    % plot lmod_st + bus P for loads listed in lmod_con...
    figure
    plot(t, bus(find(bus(:,1)==lmod_con(:,2)),6)+lmod_st)
    title([{'Modified real power loads'}; {'bus(find(bus(:,1)==lmod\_con(:,2)),6)+lmod\_st'}])
    xlabel('Time [sec]')
    ylabel('MW [PU]')
end

