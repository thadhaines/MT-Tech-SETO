% 2-machine PST simulation file for EE 5550.  Hydro speed governor on gen
% 1.  Both generators have a PSS unit.
% D. Trudnowski, Sp. 2008

bus = [ ...
% num volt angle p_gen q_gen p_load q_load G_shunt B_shunt type q_max q_min v_rated v_max v_min
  1   1.0  0.0   0.125 0     0      0      0       0       1    100   -100  20      1.5   0.5;
  2   1.0  0.0   0     0     0      0      0       0       3    0     0     20      1.5   0.5;
  3   1.0  0.0   0     0     0      0      0       0       3    0     0     20      1.5   0.5;
  4   1.0  0.0   0.25  0     0      0      0       0       2    100   -100  20      1.5   0.5
  5   1.0  0.0   0     0     0.375  0.2    0       0       3    0     0     20      1.5   0.5];

line = [ ...
% bus bus r    x    y    tapratio tapphase tapmax tapmin tapsize
  4   3   0.0  0.1  0.0  1.0      0        0      0      0; %transformer
  1   5   0.0  1.5  0.0  1.0      0        0      0      0; %line
  5   3   0.0  1.5  0.0  1.0      0        0      0      0; %line
  2   3   0.0  1.0  0.0  1.0      0        0      0      0; %line
  1   2   0.0  5.0  0.0  1.0      0        0      0      0]; %line

%Both generator parameters are from the example machine in chap 4 of
%Anderson and Fouad with a salient pole assumption.
%This model is a sub-transient model.
%From eqn (4.289), T"_qo=T'_qo if it is a 2-axis transient model.

mac_con = [ ...
% num bus base  x_l  r_a    x_d x'_d   x"_d  T'_do T"_do x_q   x'_q  x"_q  T'_qo T"_qo H      d_0  
   1   1  100   0.17 0.0    1.2 0.3    0.22  6.0   0.025 0.7   0.23  0.22  0.06  0.04  5.0    0   0   1  0.05   0.3;
   2   4  100   0.17 0.0    1.2 0.3    0.22  6.0   0.025 0.7   0.23  0.22  0.06  0.04  5.0    0   0   4  0.05   0.3];

% mac_con = [ ...
% % num bus base  x_l  r_a    x_d x'_d   x"_d  T'_do T"_do x_q   x'_q  x"_q  T'_qo T"_qo H      d_0  
%    1   1  100   0.15 0.0011 1.7 0.245  0.185 5.9   0.03  1.64  1.64  0.185 0.082 0.082 2.37   0   0   1;
%    2   4  100   0.15 0.0011 1.7 0.245  0.185 5.9   0.03  1.64  1.64  0.185 0.082 0.082 2.37   0   0   4];
  
%Exciter
%From p. 1137 of Kundur
exc_con = [...
%  type mach  T_R   K_A   T_A   T_B   T_C   Vmax  Vmin
    0    1    0     212   0     12    1     7.5   -6.0; %Fast static exciter, with TGR
    0    2    0     212   0     12    1     7.5   -6.0]; %Fast static exciter, with TGR

%PSS
%Designed using the "large-inertia, infinite-bus" method in Roger's book and
%Kundur's book.
%type gen# K  Tw T1  T2   T3  T4   max min
% pss_con = [ ...
%   1   1    80 10 0.4 0.04 0.4 0.06 0.1 -0.1;
%   1   2    80 10 0.4 0.04 0.4 0.06 0.1 -0.1];

% tg_con = [...
% %     mach wf 1/R     Tmax   Ts     Tc     T3     T4     T5
%    1   1   1  20.0    1.00   0.40   75.0   10.0  -2.4    1.2];%hydro

sw_con = [...
0        0    0    0    0    0    1/600;%sets intitial time step
15.0      3    2    0    0    0    1/600; %apply fault (fault on line 3-2)
15+2/60   0    0    0    0    0    1/600; %clear near end of fault
15+2.1/60 0    0    0    0    0    1/600; %clear far end of fault
16.0      0    0    0    0    0    1/120; % increase time step
16.1       0    0    0    0    0    1/120]; % end simulation