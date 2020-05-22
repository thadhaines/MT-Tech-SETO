% Matt Stajcar

% Input all bus, branch, and generator parameters.

function [Bus,Branch,Gen] = funBusBranchGenInputs()
%% Network Parameters.
% BUS

% Bus parameters contain the following information:
% bus number [n] - integer, ex = 5.
% bus network number - integer, ex = 1001.
% bus type: 1 = slack bus; 2 = gen bus; 3 = load bus.
% bus voltage magnitude in per unit - floating point, ex = 1.0.
% bus angle in degrees - floating point, ex = 23.25.
% rated voltage - integer, units in kilovolts, ex = 100.
% maximum voltage in per unit - floating point, ex = 1.5.
% minimum voltage in per unit - floating point, ex = 0.5.
% real power generated at bus in MW - floating point, ex = 50.0.
% reactive power generated at bus in MVAR - floating point, ex = 10.0.
% real power load at bus in MW - floating point, ex = 25.0.
% reactive power load at bus in MVAR - floating point, ex = 5.0.

% The above parameters are implemented as follows:
% Bus(52).Num = 1001;
% Bus(52).Type = 3;
% Bus(52).Vmag = 1.0;
% Bus(52).AngD = 23.25;
% Bus(52).Vrated = 100;
% Bus(52).Vmax = 1.5;
% Bus(52).Vmin = 0.5;
% Bus(52).Pgen = 0;
% Bus(52).Qgen = 0;
% Bus(52).Pload = 25.0;
% Bus(52).Qload = 5.0;

Bus(1).Num = 1;
Bus(1).Type = 2;
Bus(1).Vmag = 1.1;
Bus(1).AngD = 0;
Bus(1).Vrated = 20;
Bus(1).Vmax = 1.5;
Bus(1).Vmin = 0.5;
Bus(1).Pgen = 0.6*1.6;
Bus(1).Qgen = 0;
Bus(1).Pload = 0;
Bus(1).Qload = 0;

Bus(2).Num = 2;
Bus(2).Type = 1;
Bus(2).Vmag = 1.0;
Bus(2).AngD = 0;
Bus(2).Vrated = 20;
Bus(2).Vmax = 1.5;
Bus(2).Vmin = 0.5;
Bus(2).Pgen = 0.4*1.6;
Bus(2).Qgen = 0;
Bus(2).Pload = 0;
Bus(2).Qload = 0;

Bus(3).Num = 3;
Bus(3).Type = 3;
Bus(3).Vmag = 1.0;
Bus(3).AngD = 0;
Bus(3).Vrated = 20;
Bus(3).Vmax = 1.5;
Bus(3).Vmin = 0.5;
Bus(3).Pgen = 0;
Bus(3).Qgen = 0;
Bus(3).Pload = 1.6;
Bus(3).Qload = 0;

% BRANCH

% Branch parameters contain the following information:
% branch number [n] - integer, ex = 5.
% from bus - interger of bus number [n], ex = 1001.
% to bus - interger of bus number [n], ex = 1002.
% series impedance in per unit - R + 1j*X where R is the resistance in per  
%   unit and X is the reactance in per unit, ex = 0.1 + 1j*1.
% shunt admittance in per unit - 1j*B where B is the suseptance in per 
% unit, ex = 1j*2.
% MVA rating of branch - integer, ex = 100.

% The above paramters are implemented as follows:
% Branch(5).FmBus = 1001;
% Branch(5).ToBus = 1002;
% Branch(5).Zpu = 0.1 + 1j*1;
% Branch(5).Bpu = 1j*2;
% Branch(5).Rating = 100;

Branch(1).FmBus = 1;
Branch(1).ToBus = 3;
Branch(1).Zpu = 1j*0.152;
Branch(1).Bpu = 0;
Branch(1).Rating = 20;

Branch(2).FmBus = 2;
Branch(2).ToBus = 3;
Branch(2).Zpu = 1j*0.148;
Branch(2).Bpu = 0;
Branch(2).Rating = 20;

% GENERATOR

% Generator parameters contain the following information:
% generator number [n] - integer, ex = 5.
% bus number where generator is located - integer, ex = 1001.
% real power generated in megawatts - floating point, ex = 50.0.
% reactive power generated in megavars - floating point, ex = 10.0.
% maximum real power in megawatts - floating point, ex = 100.0.
% minimum real power in megawatts - floating point, ex = 0.0.
% maximum reactive power in megavars - floating point, ex = 35.0.
% minimum reactive power in megavars - integer, ex = -35.0.
% generator stator resistance inS per unit - floating point, ex = 0.1.
% synchronous reactance in per unit - floating point, ex = 0.2.
% transient reactance in per unit - floating point, ex = 0.2.
% subtransient reactance in per unit - floating point, ex = 0.2.

% The above parameters are implemented as follows:
% Gen(1).BusNum = 1;

Gen(1).BusNum = 1;
Gen(1).Type = 'SubTrans_Genrou';
Gen(1).R = .001;
Gen(1).Xl = .15;
Gen(1).Xd = 1.7;
Gen(1).Xd_p = 0.245;
Gen(1).Xd_pp = 0.185;
Gen(1).Td0_p = 5.9;
Gen(1).Td0_pp = 0.03;
Gen(1).Xq = 1.64;
Gen(1).Xq_p = 1.64;
Gen(1).Xq_pp = 0.185;
Gen(1).Tq0_p = 0.082;
Gen(1).Tq0_pp = 0.082;
Gen(1).Se = 1;
Gen(1).D = 0.05;
Gen(1).H = 2.37;
Gen(1).Ws = 2*pi*60;
Gen(1).dt = 1/600;

Gen(2).BusNum = 2;
Gen(2).Type = 'SubTrans_Genrou';
Gen(2).R = 0.001;
Gen(2).Xl = 0.15;
Gen(2).Xd = 1.7;
Gen(2).Xd_p = 0.245;
Gen(2).Xd_pp = 0.185;
Gen(2).Td0_p = 5.9;
Gen(2).Td0_pp = 0.03;
Gen(2).Xq = 1.64;
Gen(2).Xq_p = 1.64;
Gen(2).Xq_pp = 0.185;
Gen(2).Tq0_p = 0.082;
Gen(2).Tq0_pp = 0.082;  
Gen(2).Se = 1;
Gen(2).D = 0.02;
Gen(2).H = 5;
Gen(2).Ws = 2*pi*60;
Gen(2).dt = 1/600;

end