% script to visual compare variable and fixed data....

clear; close all

% load varible data
load('pstSETOdatane_hiskens.mat')
gv = g;
clear g;
% load fixed data
load('FTSpstSETOdatane_hiskens.mat')

busV = 16;
figure
plot(gv.sys.t(1:gv.vts.dataN), abs(gv.bus.bus_v(busV,1:gv.vts.dataN)))
hold on
plot(g.sys.t,abs(g.bus.bus_v(busV,:)) , '--')
legend('VTS','FTS','location','best')
title({ ['Bus Voltage Comparison of Bus ',int2str(busV)];'NE 39 bus 10 machine trip'})
grid on

mac = 1;

figure
plot(gv.sys.t(1:gv.vts.dataN), abs(gv.mac.mac_spd(mac,1:gv.vts.dataN)))
hold on
plot(g.sys.t,abs(g.mac.mac_spd(mac,:)) , '--')
legend('VTS','FTS','location','best')
title({ ['Machine Speed Comparison of Machine ',int2str(mac)];'NE 39 bus 10 machine trip'})
grid on

% time step size comparison
vts = zeros(gv.vts.dataN,1);
for n=2:gv.vts.dataN-1
    vts(n-1)= gv.sys.t(n)-gv.sys.t(n-1);
end

fts = zeros(size(g.sys.t,2),1);
for n=2:size(g.sys.t,2)
    fts(n-1) = g.sys.t(n) - g.sys.t(n-1);
end

figure
stairs(gv.sys.t(1:gv.vts.dataN),vts)
hold on
stairs(g.sys.t, fts)
legend('Variable ts', 'Fixed ts', 'location', 'best')
title('Time Step comparison')
ylabel('Time Step Size [seconds]')
xlabel('Time [seconds]')

