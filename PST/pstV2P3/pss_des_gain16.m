function pss_des_gain16(tw,t1,t2,t3,t4,gennum)
% function pss_des_gain16(tw,t1,t2,t3,t4,gennum)
% Function used to select the PSS gain for the 16 machine system.  This is
% run after the time constants are selected using pss_des.m.  This function
% plots the rootlocus of the system for a gain varying from 3 to 10.
% INPUTS:
%      tw = washout time constant
%      t1 = first phase lead time constant
%      t2 = first phase lag time constant
%      t3 = second phase lead time constant
%      t4 = second phase lag time constant
%      gennum = generator being analyzed
load k16a_nopss.mat; %linearize system with no pss units
G = minreal(ss(a_mat,b_vr(:,gennum),c_spd(gennum,:),0));
H=tf([tw 0],[tw 1])*tf([t1 1],[t2 1])*tf([t3 1],[t4 1]);
rlocus(H*G,-[3:10])
axis([-1 0 0 4*pi])