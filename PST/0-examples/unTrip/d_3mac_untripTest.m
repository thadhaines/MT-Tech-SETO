% 3-machine System to test un-tripping
% most model parameters from pwrmod example

bus = [ ...
% num volt  angle p_gen q_gen p_load q_load G_shunt B_shunt type q_max q_min v_rated v_max v_min
  1   1.0   0.0   1.0   0.3   0      0      0       0       1    100   -100  13.8      1.5   0.5; % Gen 1 - slack
  2   1.0   0.0   0.0   0.0   0      0      0       0       3    0     0     115      1.5   0.5;
  3   1.0   0.0   0.0   0.0   1.5    0      0       0       3    0     0     115      1.5   0.5; % load bus
  4   1.0   0.0   0.0   0.0   0      0      0       0       3    0     0     115      1.5   0.5; 
  5   1.0   0.0   0.5  0.3    0      0      0       0       2    100   -100  13.8      1.5   0.5; % Gen 2
  6   1.0   0.0   0.0    0     0      0      0       0       3    0     0     115      1.5   0.5;
  7   1.0   0.0   0.5   0      0      0      0       0       2    100   -100  13.8      1.5   0.5; % Gen 3 - small P out for tripping
  ];

line = [ ...
% bus bus r         x    y    tapratio tapphase tapmax tapmin tapsize
  1   2   0.00000  0.02500 0.0  1.0      0        1.5      .8      .05; % 'transformer'
  4   5   0.00000  0.02500 0.0  1.0      0        1.5      .8      .05; % 'transformer'
  7   6   0.0      0.005   0.00 1.0      0.       1.2       0.8     0.05; % xfmr from user manual
  %7   6   0.00000  0.02500 0.0  1.0      0        1.5      .8      .05; % 'transformer'
  2   3   0.0      0.75    0.0	1.0      0        0      0      0; % line
  4   3   0.0      0.75    0.0  1.0      0        0      0      0; % line
  6   3   0.0      0.75    0.0  1.0      0        0      0      0; % line
 % 8   6   0.0  0.2  0.0  1.0      0        0      0      0; % line
  ];

%Both generator parameters are from the example machine in chap 4 of
%Anderson and Fouad with a salient pole assumption.
%This model is a sub-transient model.
%From eqn (4.289), T"_qo=T'_qo if it is a 2-axis transient model.
mac_con = [ ...
% num bus base  x_l  r_a    x_d x'_d   x"_d  T'_do T"_do x_q   x'_q  x"_q  T'_qo T"_qo H      d_0  
   1   1  500   0.15 0.0011 1.7 0.245  0.185 5.9   0.03  1.64  1.64  0.185 0.082 0.082 2.37   0   0   1;
   2   5  200   0.15 0.0011 1.7 0.245  0.185 5.9   0.03  1.64  1.64  0.185 0.082 0.082 2.37   0   0   5;
   3   7  100   0.15 0.0011 1.7 0.245  0.185 5.9   0.03  1.64  1.64  0.185 0.082 0.082 2.37   0   0   7];
  
%Exciter
%From p. 1137 of Kundur
exc_con = [...
%  type mach  T_R   K_A   T_A   T_B   T_C   Vmax  Vmin
    0    1    0     85   0.01  12    1     7.5   -6.0; %Fast static exciter, with TGR
    0    2    0     85   0.01  12    1     7.5   -6.0; %Fast static exciter, with TGR
    0    3    0     85   0.01  12    1     7.5   -6.0; %Fast static exciter, with TGR
    ]; 

% %PSS
% %Designed using the "large-inertia, infinite-bus" method in Roger's book and
% %Kundur's book.
% %type gen# K  Tw T1  T2   T3  T4   max min
pss_con = [ ...
  1   1    80 10 0.4 0.04 0.4 0.06 0.1 -0.1;
  1   2    80 10 0.4 0.04 0.4 0.06 0.1 -0.1;
  1   3    80 10 0.4 0.04 0.4 0.06 0.1 -0.1;
  ];

tg_con = [...
%  mach wf 1/R         Tmax   Ts   Tc   T3     T4     T5
1    1  1.0 25.0    1.0     0.1     0.5 0.0     1.25  5.0 ;%user manual
1    2   1  20         1.0   0.10  5   3.0    0      0.01;%Steam
1    3   1  20         1.0   0.10  10   3.0    0      0.01;%Steam
];

 %% Load Modulation
% 
% % load_con format
% % Defines `non-conforming` loads
% % col   1 bus number
% % col   2 proportion of constant active power load
% % col   3 proportion of constant reactive power load
% % col   4 proportion of constant active current load
% % col   5 proportion of constant reactive current load
% load_con = [...
%     %   1   2       3       4       5
%         3   0       0       0       0;
%         8   0       0       0       0; % SVC
%         ]; %constant impedance
% 
% % lmod_con format
% % sets up model for load modulation
% % col	variable  						unit
% % 1  	load modulation number 
% % 2  	bus number 
% % 3  	modulation base MVA  			MVA
% % 4  	maximum conductance lmod_max 	pu
% % 5  	minimum conductance lmod_min 	pu
% % 6  	regulator gain K  				pu
% % 7  	regulator time constant T_R  	sec
% lmod_con = [ ...
%   % 1   2   3       4       5       6       7
%     1   3   100     1.0     0.0     1.0     0.1;];
% %svc
% % col 1           svc number
% % col 2           bus number
% % col 3           svc base MVA
% % col 4           maximum susceptance Bcvmax(pu)
% % col 5           minimum susceptance Bcvmin(pu)
% % col 6           regulator gain
% % col 7		  regulator time constant (s)
% 
% svc_con = [1  8  100  1  -1  10  0.05];
% 

%% Switching
%Switching file defines the simulation control
% row 1 col1  simulation start time (s) (cols 2 to 6 zeros)
%     col7  initial time step (s)
% row 2 col1  fault application time (s)
%     col2  bus number at which fault is applied
%     col3  bus number defining far end of faulted line
%     col4  zero sequence impedance in pu on system base
%     col5  negative sequence impedance in pu on system base
%     col6  type of fault  - 0 three phase
%                  - 1 line to ground
%                  - 2 line-to-line to ground
%                  - 3 line-to-line
%                  - 4 loss of line with no fault
%                  - 5 loss of load at bus
%     col7  time step for fault period (s)
% row 3 col1  near end fault clearing time (s) (cols 2 to 6 zeros)
%     col7  time step for second part of fault (s)
% row 4 col1  far end fault clearing time (s) (cols 2 to 6 zeros)
%     col7  time step for fault cleared simulation (s)
% row 5 col1  time to change step length (s)
%     col7  time step (s)
sw_con = [...
0        0    0    0    0    0    1/120; % sets intitial time step
2.0      6    1    0    0    6    1/120; % no switching action - trip gen via mac_trip_logic
15.0     0    0    0    0    0    1/120; % reconnect gen
20.0     0    0    0    0    0    1/120; % SS reached, un-trip gen via mac_trip_logic
25.0      0    0    0    0    0    1/120; % add Gov droop control
35.0      0    0    0    0    0    1/120; % add exciter
45.0      0    0    0    0    0    1/120; % ramp tg sig to gov Pref
65.0      0    0    0    0    0    1/120; % set tg sig to Pref
70.0      0    0    0    0    0    1/120; % set Pref in gov, remove sig
80.0      0    0    0    0    0    1/120; % start ramping exciter in
120.0     0    0    0    0    0    1/120; % clear far end of fault - end simulation
    ]; % end simulation

%% solver_con format
% A cell with a solver method in each row corresponding to the specified
% 'time blocks' defined in sw_con
%
% Valid solver names:
% huens - Fixed time step default to PST
% ode113 - works well during transients, consistent # of slns, time step stays relatively small
% ode15s - large number of slns during init, time step increases to reasonable size
% ode23 - realtively consisten # of required slns, timstep doesn't get very large
% ode23s - many iterations per step - not efficient...
% ode23t - occasionally hundereds of iterations, sometimes not... decent
% ode23tb - similar to 23t, sometimes more large sln counts
% 
% solver_con ={ ...
%     'ode113'; % pre fault - fault
%     'ode23'; % fault - post fault 1
%     'huens'; % post fault 1 - post fault 2
%     'ode23';
%     'ode23';
%     };
