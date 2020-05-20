%% Ybus Tester.

% Matt Stajcar

clear all; close all; clc

%% Run PST load flow solver.

% d1m_classical
lfdemo

pf.bus = bus_sol(:,1);
pf.Vmag = bus_sol(:,2);
pf.VangD = bus_sol(:,3);
pf.Pg = bus_sol(:,4);
pf.Qg = bus_sol(:,5);
pf.Pl = bus_sol(:,6);
pf.Ql = bus_sol(:,7);

pf.VangR = pf.VangD * pi/180;
pf.Vcart = pf.Vmag.*(cos(pf.VangR) + 1j*sin(pf.VangR));
pf.Sg = pf.Pg + 1j*pf.Qg;
pf.Sl = pf.Pl + 1j*pf.Ql;
pf.S = pf.Sg + pf.Sl;

pf.Ig = conj(pf.S./pf.Vcart);
Iest = Ybus.Y*pf.Vcart;
Scheck = pf.Vcart.*conj(pf.Ig);
    



