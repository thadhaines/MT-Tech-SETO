% 1-machine infinite-bus system for EE 5550.  PST simulation file
% Dan Trudnowski, 2010

bus = [ ...
% num volt angle p_gen q_gen p_load q_load G_shunt B_shunt type q_max q_min v_rated v_max v_min
  1   1.0  0.0   -0.25 0     0      0      0       0       1    100   -100  20      1.5   0.5;
  2   1.0  0.0   0     0     0      0      0       0       3    0     0     20      1.5   0.5;
  3   1.0  0.0   0     0     0      0      0       0       3    0     0     20      1.5   0.5;
  4   1.0  0.0   0.25  0     0      0      0       0       2    100   -100  20      1.5   0.5];

line = [ ...
% bus bus r    x    y    tapratio tapphase tapmax tapmin tapsize
  4   3   0.0  0.1  0.0  1.0      0        0      0      0; %transformer
  1   3   0.0  3.0  0.0  1.0      0        0      0      0; %line
  2   3   0.0  1.0  0.0  1.0      0        0      0      0; %line
  1   2   0.0  5.0  0.0  1.0      0        0      0      0]; %line

%Gen 2 is a classical machine from chap. 4 of Anderson and Fouad
mac_con = [ ...
% num bus base  x_l  r_a   x_d x'_d   x"_d  T'_do T"_do x_q   x'_q  x"_q  T'_qo T"_qo H      d_0  d_1  bus
   1   1  100   0    0     0   0.0003 0     0     0     0     0     0     0     0     2500   0     0   1; %Infinite bus
   2   4  100   0    0     0   0.245  0     0     0     0     0     0     0     0     2.37   0     0   4];

%Switching file defines the simulation control
% row 1 
%   col1 simulation start time (s) (cols 2 to 6 = zeros), 
%   col7 initial time step (s)
% row 2 
%   col1 fault application time (s)
%   col2 bus number at which fault is applied
%   col3 bus number defining far end of faulted line
%   col4 zero sequence impedance in pu on system base
%   col5 negative sequence impedance in pu on system base
%   col6 type of fault 
%   - 0 three phase
%   - 1 line to ground
%   - 2 line-to-line to ground
%   - 3 line-to-line
%   - 4 loss of line with no fault
%   - 5 loss of load at bus
%   col7 time step for fault period (s)
% row 3 
%   col1 near end fault clearing time (s) (cols 2 to 6 zeros)
%   col7 time step for second part of fault (s)
% row 4 
%   col1 far end fault clearing time (s) (cols 2 to 6 zeros)
%   col7 time step for fault cleared simulation (s)
% row 5 
%   col1 time to change step length (s)
%   col7 time step (s)  
sw_con = [...
0        0    0    0    0    0    1/600;%sets intitial time step
1.0      3    2    0    0    0    1/600; %apply fault (fault on line 3-2)
1+2/60   0    0    0    0    0    1/600; %clear near end of fault
1+2.1/60 0    0    0    0    0    1/600; %clear far end of fault
2.0      0    0    0    0    0    1/120; % increase time step
20       0    0    0    0    0    1/120]; % end simulation