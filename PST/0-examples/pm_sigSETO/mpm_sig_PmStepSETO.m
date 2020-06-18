function f = mpm_sig(t, k)
% Syntax: f = mpm_sig(t,k)
% 1:19 PM 15/08/97
% defines modulation signal for generator mechanical power
global pm_sig n_pm
global g
f = 0;
if n_pm~=0
%   pm_sig(:,k) = zeros(n_pm,1);
   if g.sys.t(k)> 1.0
      pm_sig(1,k) = -0.025;
   end
   %  pm_sig(:,k) = zeros(n_pm,1);
   %pm_sig(1,k) = 0.01;
   %pm_sig(2,k) = -0.01;
  %end
end
return
