function mpm_sig( k)
% Syntax: f = mpm_sig(t,k)
% 1:19 PM 15/08/97
% defines modulation signal for generator mechanical power
%global pm_sig n_pm
global g

if g.mac.n_pm~=0
    if g.sys.t(k)> 1.0
        g.mac.pm_sig(1,k) = -0.025;
    end
    if g.sys.t(k)> 5
        g.mac.pm_sig(1,k) = 0;
    end
end
return
