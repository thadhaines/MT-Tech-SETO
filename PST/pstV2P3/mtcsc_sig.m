function f = mtcsc_sig(t,k)
% Syntax: f = mtcsc_sig(t,k)
% 4:39 PM 15/08/97
% defines modulation signal for tcsc control
global tcsc_sig n_tcsc mac_spd
f=0; %dummy variable
if n_tcsc ~=0
  tcsc_sig(:,k) = zeros(n_tcsc,1);
end