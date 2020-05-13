function f = mpm_sig(t,k)
% Syntax: f = mpm_sig(t,k)
% 1:19 PM 15/08/97
% defines modulation signal for generator mechanical power
global pm_sig n_pm
f=0; %dummy variable
if n_pm~=0
   pm_sig(:,k) = zeros(n_pm,1);
%    if t>=1 & t<2.5
%        pm_sig(7,k) = 0.0; %Pulse generator 7
%    end
end