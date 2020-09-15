function mpm_sig(t, k)
% Syntax: mpm_sig(t,k)
% 1:19 PM 15/08/97
% defines modulation signal for generator mechanical power

global n_pm pm_sig

if n_pm~=0
    if t > 1.0
        pm_sig(1,k) = -0.025;
    end
    if t > 5
        pm_sig(1,k) = 0;
    end
end
return
