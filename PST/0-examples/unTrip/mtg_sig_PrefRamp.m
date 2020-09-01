function mtg_sig(k)
% MTG_SIG Defines modulation signal for turbine power reference
% Syntax: mtg_sig(k)
%
global g
% actions to return a generator back on line

%% ramp Pref near to original value
if abs(g.sys.t(k)-45) < 1e-6
    disp(['MTG_SIG:  ramping gov Pref via tg_sig at t = ', num2str(g.sys.t(k))])
end
if g.sys.t(k)>= 45 && g.sys.t(k)< 65 %
    g.tg.tg_sig(3,k) = (g.sys.t(k)-45)*(0.5003)/20; % 25 second ramp up
end

%% set signal near to pref,
if abs(g.sys.t(k)-65) < 1e-6
    disp(['MTG_SIG:  sig ramp done, setting sig at t = ', num2str(g.sys.t(k))])
end
if g.sys.t(k)>= 65 && g.sys.t(k)<70
    g.tg.tg_sig(3,k) = (0.5003);
end

%% set pref, remove signal
if abs(g.sys.t(k)-70) < 1e-6
    g.tg.tg_sig(3,k) = 0; % remove Pref sig
    g.tg.tg_pot(3,5) = 0.5003; % set Pref
    reInitGov(3,k)
    disp(['MTG_SIG:  setting Pref, removing sig at t = ', num2str(g.sys.t(k))])
end
end% end function
