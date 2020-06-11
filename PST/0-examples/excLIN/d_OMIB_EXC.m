%%  d_OneMacInfBus.m
%   Thad Haines         EELE 5550
%   Program Purpose:    Create bus, line, mac_con and sw_con for sub transient PST Model. 
%                       DC exciter, steam gov, all bus same V base, lossless lines, 

%   History:
%   03/20/18    21:12   init
%   03/20/18    21:19   copied from base, added a DC exciter (type 1)
%   05/19/20    10:26   Modified for use in exploring PST operation
%   05/21/20    08:24   Added optional Fbase and Sbase variables
%   06/11/20    09:22   Removed load modulation for linear comparison

%% Optional System data, else assumed 60 Hz, 100 MVA
%Fbase = 50; % Hz
%Sbase = 50; % MW

%% bus data 
% bus array format:
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
%     bus_type  - 1, swing bus
%               - 2, generator bus (PV bus)
%               - 3, load bus (PQ bus)
% col11 q_gen_max(pu)
% col12 q_gen_min(pu)
% col13 v_rated (kV)
% col14 v_max pu
% col15 v_min pu
bus = [ ...
%   1   2   3       4       5       6       7       8   9   A   B       C       D        E   F
%busNum V delta     P_gen   Q_gen P_load   Q_load   G   B  type Q_max  Q_minVrated_(kv) Vmax Vmin
    1   1   0       -.5     0       0       0       0   0   1   100     -100  	20       2   0.5; %slack / infinite bus
    2   1   0       0       0       0.3     0       0   0   3   0       0       20       2   0.5;
    3   1   0       0       0       0       0       0   0   3   0       0       20       2   0.5;
    4   1   0       0.25    0       0       0       0   0   2   100     -100  	20       2   0.5];    %gen
   
%% line data
% line array format:
% col1 from bus
% col2 to bus
% col3 resistance(pu)
% col4 reactance(pu)
% col5 line charging(pu)
% col6 tap ratio
% col7 tap phase
% col8 tapmax
% col9 tapmin
% col10 tapsize

line = [ ...
%   1   2   3   4       5   6   7       8   9   A   
%From   To  R   X       LC tp/ Phase   max min size  
    1   2   0   5       0   1   0       0   0   0;
    1   3   0   3       0   1   0       0   0   0;
    2   3   0   1       0   1   0       0   0   0;
    3   4   0   0.1     0   1   0       0   0   0]; % 'transformer'?
    
%% Machine data 
% mac_con array format:
% 1. machine number,
% 2. bus number,
% 3. base mva,
% 4. leakage reactance x_l(pu),
% 5. resistance r_a(pu),
% 6. d-axis sychronous reactance x_d(pu),
% 7. d-axis transient reactance x'_d(pu),
% 8. d-axis subtransient reactance x"_d(pu),
% 9. d-axis open-circuit time constant T'_do(sec),
% 10. d-axis open-circuit subtransient time constant T"_do(sec),
% 11. q-axis sychronous reactance x_q(pu),
% 12. q-axis transient reactance x'_q(pu),
% 13. q-axis subtransient reactance x"_q(pu),
% 14. q-axis open-circuit time constant T'_qo(sec),
% 15. q-axis open circuit subtransient time constant T"_qo(sec)
% 16. inertia constant H(sec),
% 17. damping coefficient d_o(pu),
% 18. dampling coefficient d_1(pu),
% 19. bus number
%

mac_con = [ ...
%   1   2   3 [MVA] 4       5       6       7       8       9       A       B       C       D       E       F       16      17      18      19
% Mac_# b_# Sbase   x_l     r_a     x_d     x'_d    x"_d    T'_do   T"_do   x_q     x'_q    x"_q	T'_qo   T"_qo   H       d_o     d_1     b_#
    1   1   100     0.15    0.0011  1.7     0.245   0.185   5.9     0.03    1.64    1.64    0.185   0.082   0.082   5.0    0       0       1;  % G1 - double size
    2   4   100     0.15    0.0011  1.7     0.245   0.185   5.9     0.03    1.64    1.64    0.185   0.082   0.082   2.37    0       0       4];   % G2
                                            % infinite bus really small x'_d
                                            % if x'_d = zero -> crash?
                                            
%% Exciter data
% exc_con format:
% column    data  									unit
% 1  		exciter type  							3 for ST3
% 2  		machine number 
% 3         input filter time constant T_R			sec
% 4         voltage regulator gain K_A 
% 5  		voltage regulator time constant T_A		sec
% 6  		voltage regulator time constant T_B		sec
% 7  		voltage regulator time constant T_C		sec
% 8  		maximum voltage regulator output V_Rmax pu
% 9  		minimum voltage regulator output V_Rmin	pu
% 10		maximum internal signal V_Imax 			pu
% 11		minimum internal signal V_Imin			pu
% 12		first state regulator gain K_J
% 13        potential circuit gain coefficient K_P
% 14        potential circuit phase angle qP 		degrees
% 15        current circuit gain coefficient K_I
% 16		potential source reactance X_L 			pu
% 17		rectifier loading factor K_C 
% 18		maximum field voltage E_fdmax 			pu
% 19		inner loop feedback constant K_G 
% 20        max inner loop voltage feedback V_Gmax  pu

exc_con = [ ...
%   1       2       3       4       5       6       7       8       9      
%   Type    mach #  T_R     K_A     T_A     T_B     T_C     V_Rmax  V_Rmin  
    0       1       0       100      0.01    12.0    1.0     7.5     -6;   
%    0       2       0       150      0.01    12.0    1.0     7.5     -6; 
];

%% Governor data
% tg_con matrix format
%column        data         unit
%  1    turbine model number (=1)   
%  2    machine number  
%  3    speed set point   wf        pu
%  4    steady state gain 1/R       pu
%  5    maximum power order  Tmax   pu on generator base
%  6    servo time constant   Ts    sec
%  7    governor time constant  Tc  sec
%  8    transient gain time constant T3 sec
%  9    HP section time constant   T4   sec
% 10    reheater time constant    T5    se

tg_con = [ ...
%   1   2   3   4       5       6       7       8       9       A
%               1/R     Tmax    Ts      Tc      T3      T4      T5
    1   1   1   1/0.05  1.00    0.4     45.00   5.00    -1.00   0.5 % hydro Gov
    1   2   1   1/0.05  1.00    0.04    0.20    0.0     1.50    5.00 %steam Gov
];

%% Switching file 
% row 1 col1 simulation start time (s) (cols 2 to 6 zeros)
%       col7 initial time step (s)
% row 2 col1 fault application time (s)
%       col2 bus number at which fault is applied
%       col3 bus number defining far end of faulted line
%       col4 zero sequence impedance in pu on system base
%       col5 negative sequence impedance in pu on system base
%       col6 type of fault - 0 three phase
%           - 1 line to ground
%           - 2 line-to-line to ground
%           - 3 line-to-line
%           - 4 loss of line with no fault
%           - 5 loss of load at bus
%           - 6 no action
%       col7 time step for fault period (s)
% row 3 col1 near end fault clearing time (s) (cols 2 to 6 zeros)
%       col7 time step for second part of fault (s)
% row 4 col1 far end fault clearing time (s) (cols 2 to 6 zeros)
%       col7 time step for fault cleared simulation (s)
% row 5 col1 time to change step length (s)
%       col7 time step (s)
%
% row n col1 finishing time (s) (n indicates that intermediate rows may be inserted)

sw_con = [ ...
%   1                   2   3   4   5   6   7   
    0                   0   0   0   0   0   0.01;     % time start, ts set
    1                   3   2   0   0   6   0.01;  % No Action
    1+(2/60)            0   0   0   0   0   0.01;  % 
    1+(2/60)+(1/600)    0   0   0   0   0   0.01;   % 
    20                   0   0   0   0   0   0.01;   % end
];    
  
%% Load Modulation

% load_con format
% Defines `non-conforming` loads
% col   1 bus number
% col   2 proportion of constant active power load
% col   3 proportion of constant reactive power load
% col   4 proportion of constant active current load
% col   5 proportion of constant reactive current load
%load_con = [...
    %   1   2       3       4       5
%        2   1.0       1.0       0       0;]; %constant power
%
% lmod_con format
% sets up model for load modulation
% col	variable  						unit
% 1  	load modulation number 
% 2  	bus number 
% 3  	modulation base MVA  			MVA
% 4  	maximum conductance lmod_max 	pu
% 5  	minimum conductance lmod_min 	pu
% 6  	regulator gain K  				pu
% 7  	regulator time constant T_R  	sec
% lmod_con = [ ...
%   % 1   2   3       4       5       6       7
%     1   2   100     1.0     0.0     1.0     0.1;];

