% Governor Runner.

clear
close all
clc

% Gov params.
R = 0.05;
T1 = 0.5;
T2 = 3.0;
T3 = 10;
Dt = 0;

tstart = 0;
dt = .001;
tstop = 600-dt;

t = (tstart:dt:tstop)';

% Pref = [t,.25*ones(length(t),1)];
deltaW = [t,zeros(length(t),1)];

Wchange = 0.01*ones(400e3,1);
deltaW(200e3:600e3-1,2) = Wchange;

sim('GovModel')
figure
plot(t,Pmech)

sim('GovModel02')
figure
plot(t,Pmech,'r')


