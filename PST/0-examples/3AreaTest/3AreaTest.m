% Three Area test - with load step in area 1 capabilities
%

enableGov = true;
enableAGC = false;
conditionalAGC = 0; % 1 or 0

disp('Three areas, each with one gen bus and one load bus')

%% bus data format
%{
% bus: 
% col1 number
% col2 voltage magnitude(pu)
% col3 voltage angle(degree)
% col4 p_gen(pu)
% col5 q_gen(pu),
% col6 p_load(pu)
% col7 q_load(pu)
% col8 G shunt(pu)
% col9 B shunt(pu)
% col10 bus_type
%       bus_type - 1, swing bus
%               - 2, generator bus (PV bus)
%               - 3, load bus (PQ bus)
% col11 q_gen_max(pu)
% col12 q_gen_min(pu)
% col13 v_rated (kV)
% col14 v_max  pu
% col15 v_min  pu
%}
%           v       ang   pgen   qgen  pL    qL    G     B    type Qmx  Qmn   kv   Vmx  Vmn
bus = [ 11  1.00    0.0   0.09   0.00  0.00  0.00  0.00  0.00 1  99.0 -99.0 115.0   1.1  .9;
        12  1.00    0.0   0.00   0.00  0.08  0.00  0.00  0.00 3  0.0   0.0  115.0  1.05 .95;
        % middle 
        21  1.00    0.0   0.04   0.00  0.00  0.00  0.00  0.00 2  10.0  -5.0  115.0   1.1  .9;
        22  1.00    0.0   0.00   0.00  0.05  0.00  0.00  0.00 3  0.0   0.0  115.0  1.05 .95;
        % south
        31  1.00    0.0   0.17   0.00  0.00  0.00  0.00  0.00 2  10.0  -5.0  115.0   1.1  .9;
        32  1.00    0.0   0.00   0.00  0.17  0.00  0.00  0.00 3  0.0   0.0  115.0  1.05 .95 ];
   
%% area_def data format
% should contain same number of rows as bus array (i.e. all bus areas defined)
% col 1 bus number
% col 2 area number
area_def = [ ...
            11  1;
            12  1;
            21  2;
            22  2;
            31  3;
            32  3];
   
%% line data format
% line: from bus, to bus, resistance(pu), reactance(pu),
%       line charging(pu), tap ratio, tap phase, tapmax, tapmin, tapsize

line = [...
11  12  0.0     0.005    0.00    1.0  0. 0.  0.  0.;
21  22  0.0     0.005    0.00    1.0  0. 0.  0.  0.;
31  32  0.0     0.005    0.00    1.0  0. 0.  0.  0.;

12  21  0.0     0.0167   0.00    1.0  0. 0.  0.  0.;
22  31  0.0     0.0167   0.00    1.0  0. 0.  0.  0.;
32  11  0.0     0.0167   0.00    1.0  0. 0.  0.  0.];

%% Line Monitoring
% Each value corresponds to an array index in the line array.
% Complex current and power flow on the line will be calculated and logged during simulation

lmon_con = [4, 5, 6];

%% Machine data format
%{
%       1. machine number,
%       2. bus number,
%       3. base mva,
%       4. leakage reactance x_l(pu),
%       5. resistance r_a(pu),
%       6. d-axis sychronous reactance x_d(pu),
%       7. d-axis transient reactance x'_d(pu),
%       8. d-axis subtransient reactance x"_d(pu),
%       9. d-axis open-circuit time constant T'_do(sec),
%      10. d-axis open-circuit subtransient time constant T"_do(sec),
%      11. q-axis sychronous reactance x_q(pu),
%      12. q-axis transient reactance x'_q(pu),
%      13. q-axis subtransient reactance x"_q(pu),
%      14. q-axis open-circuit time constant T'_qo(sec),
%      15. q-axis open circuit subtransient time constant T"_qo(sec),
%      16. inertia constant H(sec),
%      17. damping coefficient d_o(pu),
%      18. dampling coefficient d_1(pu),
%      19. bus number
%
% note: all the following machines use sub-transient model
%}
mac_con = [ ...

1 11 900 0.200  0.0025  1.8  0.30  0.25 8.00  0.03...
                       1.7  0.55  0.25 0.4   0.05...
  6.5  0  0  1;
2 21 900 0.200  0.0025  1.8  0.30  0.25 8.00  0.03...
                       1.7  0.55  0.25 0.4   0.05...
  6.5  0  0  2;
3 31 900 0.200  0.0025  1.8  0.30  0.25 8.00  0.03...
                       1.7  0.55  0.25 0.4   0.05...
  6.5  0  0  11];


%% governor model
%{
% tg_con matrix format
%column	       data			unit
%  1	turbine model number (=1)	
%  2	machine number	
%  3	speed set point   wf		pu
%  4	steady state gain 1/R		pu
%  5	maximum power order  Tmax	pu on generator base
%  6	servo time constant   Ts	sec
%  7	governor time constant  Tc	sec
%  8	transient gain time constant T3	sec
%  9	HP section time constant   T4	sec
% 10	reheater time constant    T5	sec
%}
if enableGov
    tg_con = [...
    1  1  1  10.0  1.0  0.1  0.5 0.0 1.25 5.0;
    1  2  1  10.0  1.0  0.1  0.5 0.0 1.25 5.0;
    1  3  1  10.0  1.0  0.1  0.5 0.0 1.25 5.0];
end

%% non-conforming load
%{
% col 1           bus number
% col 2           fraction const active power load
% col 3           fraction const reactive power load
% col 4           fraction const active current load
% col 5           fraction const reactive current load
%}
load_con = [...
12   0  0   0  0;
22   0  0   0  0;
32   0  0   0  0;
]; 

%% lmod_con format
%{
% sets up model for load modulation
% col	variable  						unit
% 1  	load modulation number 
% 2  	bus number 
% 3  	modulation base MVA  			MVA
% 4  	maximum conductance lmod_max 	pu
% 5  	minimum conductance lmod_min 	pu
% 6  	regulator gain K  				pu
% 7  	regulator time constant T_R  	sec
%}
lmod_con = [ ...
  % 1   2   3       4       5       6       7
    1   12   100     10.0     0.0     1.0     0.004;
   % 1   14   100    1.0     0.0     1.0     0.004 ;
   ];

%% AGC definition
%{ 
Experimental model definition akin to Trudnowski experimental code.
Each agc(x) has following fields:
    area        - Area number / controlled area
    startTime   - Time of first AGC signal to send
    actionTime  - Interval of time between all following AGC signals
    gain        - Gain of output ACE signal
    Btype       - Fixed frequency bias type (abs, percent of max capacity...)
        0 - absolute - Input B value is set as Frequency bias (positive MW/0.1Hz)
        1 - percent of max area capacity
    B           - Fixed frequency bias Value
    Kbv         - Varaible frequency bias gain used to gain B as B(1+kBv*abs(delta_w))
    condAce     - Conditional ACE flag
        0 - Conditional ACE not considered
        1 - ace2dist updated only if sign matches delta_w (i.e. in area event)

    (PI Filter Values)
    Kp          - Proportional gain
    a           - ratio between integral and proportional gain (placement of zero)

    ctrlGen_con - Controlled generator information (see format note below)
%}
agc(1).area = 1;
agc(1).startTime = 25;
agc(1).actionTime = 15;
agc(1).gain = 2; % gain of output signal
agc(1).Btype = 1; % per max area capacity
agc(1).B = 1;
agc(1).Kbv = 0; % no variable bias
agc(1).condAce = conditionalAGC; % conditional ACE
agc(1).Kp = 0.04;
agc(1).a = 0.001;
agc(1).ctrlGen_con = [ ...
    % ctrlGen_con Format:
    %col 1 gen bus
    %col 2 participation Factor
    %col 3 low pass filter time constant [seconds]
    11, 1, 2;
    ];

agc(2)=agc(1); % duplicate most settings from AGC 1 to AGC 2
agc(2).area = 2;
agc(2).ctrlGen_con = [...
%    col 1 gen bus
%    col 2 participation Factor
    %col 3 low pass filter time constant [seconds]
    21, 1, 2;
    ];

agc(3)=agc(1); % duplicate most settings from AGC 1 to AGC 3
agc(3).area = 3;
agc(3).ctrlGen_con = [...
%    col 1 gen bus
%    col 2 participation Factor
    %col 3 low pass filter time constant [seconds]
    31, 1, 2;
    ];

if ~enableAGC
    % set gains to zero
    agc(1).gain = 0;
    agc(2).gain = 0;
    agc(3).gain = 0;
end

%% Switching file defines the simulation control
%{
% row 1 col1  simulation start time (s) (cols 2 to 6 zeros)
%       col7  initial time step (s)
% row 2 col1  fault application time (s)
%       col2  bus number at which fault is applied
%       col3  bus number defining far end of faulted line
%       col4  zero sequence impedance in pu on system base
%       col5  negative sequence impedance in pu on system base
%       col6  type of fault  - 0 three phase
%                            - 1 line to ground
%                            - 2 line-to-line to ground
%                            - 3 line-to-line
%                            - 4 loss of line with no fault
%                            - 5 loss of load at bus
%       col7  time step for fault period (s)
% row 3 col1  near end fault clearing time (s) (cols 2 to 6 zeros)
%       col7  time step for second part of fault (s)
% row 4 col1  far end fault clearing time (s) (cols 2 to 6 zeros)
%       col7  time step for fault cleared simulation (s)
% row 5 col1  time to change step length (s)
%       col7  time step (s)
% row n col1 finishing time (s)  (n indicates that intermediate rows may be inserted)
%}

ts = 0.004;
sw_con = [...
0    0    0    0    0    0    ts;   % sets intitial time step
0.5  12  21    0    0    6    ts;   % Do Nothing
0.75  0  0      0    0    0    ts;   % Do Nothing
30.0  0      0    0    0    0    ts;   % Do Nothing
];   % end simulation

clear ts enableExciters enableGov enablePSS enableSVC conditionalAGC