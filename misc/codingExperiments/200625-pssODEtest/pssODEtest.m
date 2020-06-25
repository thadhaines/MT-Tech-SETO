%% test to use ode solver to step PST model
close all
clear
clc
%% globals that are pss related
global  pss_con pss_pot mac_con mac_int exc_con; 
global pss_idx n_pss pss_sp_idx pss_p_idx pss_mb_idx pss_exc_idx;
global pss_T  pss_T2 pss_T4 pss_T4_idx  pss_noT4_idx;

%% globals that are user supplied or relate to a machine...
%           1   2   3   4       5   6   7       8       9       10
pss_con = [ 1  	1   20  2     0.25 0.04 0.2   0.03      1.0     -1.0];
exc_con = [...
%   1   2   3     4     5     6     7    8     9     10    11    12  13  14  15  16  17 18?    19  20   
%  type num T_R   K_A   T_A   T_B   T_C  Vrmax Vrmin VImax VImin KJ  KP  qP  KI  XL  KC Efdmax KG  VGmax
    3    1  0.0   200   0.02  10.0  1.0  5.0   -5.0  0.1   -0.1  1.0 1   0   0   0   0  5.0    0  5.0];
mac_int = 1;
d_0S = 0; %Salient damping
d_0R = 0; %Round rotor damping
d_1 = 1; %Pmech damping (Not used in new mac_sub)
mac_con = [
% 1   2   3      4x    5y   6y   7y    8y    9y     10y    11y   12y   13y   14y   15y    16y   17    18   19  20x    21x
% num bus base   x_l   r_a  x_d  x'_d  x"_d  T'_do  T"_do  x_q   x'_q  x"_q  T'_qo T"_qo  H     d_0   d_1  bus s(1.0) s(1.2)   
   1   1  4500   0.17  0.0  1.2  0.3   0.22  6.0    0.025  0.7   0.23  0.22  0.06  0.04   5.0   d_0S  d_1   1  0.05   0.3; %Hydro, salient pole
];
%% MATLAB model
block1 = tf([pss_con(3)*pss_con(4), 0],[pss_con(4), 1])
block2 = tf([pss_con(5), 1],[pss_con(6), 1])
block3 = tf([pss_con(7), 1],[pss_con(8), 1])

G= block1*block2*block3;
tL = 0:.001:5;
modSig = zeros(size(tL,1),1);
modSig(tL>=1) = .001; % very small input to avoid limiter
yL = lsim(G,modSig,tL);
figure
plot(tL,yL)

%% PST - INCOMPLETE
% used to make indexes of pss
pss_indx

%% initialize zeros of pss variables...
k = size(tL,2);
z_pss = zeros(1,k);
if n_pss~=0
    z_pss = zeros(n_pss,k);
end

pss1 = z_pss; 
pss2 = z_pss; 
pss3 = z_pss;
dpss1 = z_pss; 
dpss2 = z_pss; 
dpss3 = z_pss;

%%
bus= 'dummy';
pss23(0,1,bus,0)
