function mpm_sig( k)
% Syntax: f = mpm_sig(t,k)

% defines modulation signal for generator mechanical power

% NOTE: pmech is altered directly, though pm_sig could also be used / was 'the way'
global g

if g.mac.n_pm~=0
    if g.sys.t(k)>= 40 && g.sys.t(k)< 65 %
        g.mac.pm_sig(3,k) = (g.sys.t(k)-40)*0.5/25; % 25 second ramp up
        %g.mac.pmech(3,k) = (g.sys.t(k)-40)*0.5/25; % 25 second ramp up
    end
    if g.sys.t(k)>= 65
        g.mac.pm_sig(3,k) = .5;
        %g.mac.pmech(3,k) = .5; % setting pmech back to original setting
    end
end
return
