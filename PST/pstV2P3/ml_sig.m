function f = ml_sig(t,k)
% Syntax: f = ml_sig(t,k)
%4:40 PM 15/08/97
% defines modulation signal for lmod control
global lmod_sig n_lmod bus_v mac_spd
f=0; %dummy variable

lmod_sig(:,k) = zeros(size(lmod_sig(:,k))); %Initialize current value to zero
