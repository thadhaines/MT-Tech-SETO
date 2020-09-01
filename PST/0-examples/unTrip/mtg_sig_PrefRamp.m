function mtg_sig(k)
% MTG_SIG Defines modulation signal for turbine power reference
% Syntax: mtg_sig(k)
%
global g
% actions to return a generator back on line
% ramp pref instead of tg sig

%% ramp Pref near to original value
if abs(g.sys.t(k)-45) < 1e-6
    disp(['MTG_SIG:  ramping gov Pref at t = ', num2str(g.sys.t(k))])
end
if g.sys.t(k)>= 45 && g.sys.t(k)< 65 %
    g.tg.tg_pot(3,5) = (g.sys.t(k)-45)*(0.5003)/20; % ramp reference
end

%% set signal near to pref,
if abs(g.sys.t(k)-65) < 1e-6
    disp(['MTG_SIG:  Pref ramp done, setting Pref at t = ', num2str(g.sys.t(k))])
    g.tg.tg_pot(3,5) = (0.5003);
end

end% end function
