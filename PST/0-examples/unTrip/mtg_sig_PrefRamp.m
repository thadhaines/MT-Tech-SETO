function mtg_sig(k)
% MTG_SIG Defines modulation signal for turbine power reference
% Syntax: mtg_sig(k)
%
global g
% actions to return a generator back on line

% ramp Pref near to original value
if abs(g.sys.t(k)-45) < 1e-5
    disp('ramping gov Pref via tg_sig')
end
if g.sys.t(k)>= 45 && g.sys.t(k)< 65 %
    g.tg.tg_sig(3,k) = (g.sys.t(k)-45)*(0.5)/20; % 25 second ramp up
end

% set signal near to pref,
if abs(g.sys.t(k)-65) < 1e-5
    disp('sig ramp done, setting sig')
end
if g.sys.t(k)>= 65 && g.sys.t(k)<70
    g.tg.tg_sig(3,k) = (0.5);
end

% set pref, remove signal
if abs(g.sys.t(k)-70) < 1e-5
    g.tg.tg_sig(3,k) = 0; % remove Pref sig
    g.tg.tg_pot(3,5) = 0.5; % set Pref
    reInitGov(3,k) 
    disp('setting Pref, removing sig')
end

end
