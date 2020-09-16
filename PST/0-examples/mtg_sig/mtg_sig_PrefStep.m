function mtg_sig(t, k)
% MTG_SIG Defines modulation signal for turbine power reference
% Syntax: f = mtg_sig(t, k)

global tg_sig

if t(k) > 0.5
    tg_sig(1,k) = 0.05; % modify first load only
end
if t(k) > 5
    tg_sig(1,k) = -0.05; % to counter lmod above...
end
if t(k) > 10
    tg_sig(1,k) = 0; % to counter lmod above...
end
return
