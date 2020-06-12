function mtg_sig(k)
% MTG_SIG Defines modulation signal for turbine power reference
% Syntax: f = mtg_sig(k)
%
%   History:
%   Date        Time    Engineer        Description
%   07/xx/98    12:37   Graham Rogers   Version 1
%   06/05/20    16:21   Thad Haines     V2 - using global g, no t passed
%   in, no dummy varibale f passed out. Default behavior commented out.
%

global g

%fprintf('%4.4f \t %d\n', t(k), k); % DEBUG
% probably less sloppy way to do this...
if g.sys.t(k) > 0.5
    %g.lmod.lmod_sig(1,k) = 0.25; % modify first load only
    % Pref step
    g.tg.tg_sig(1,k) = 0.05; % to counter lmod above...
end
if g.sys.t(k) > 5
    g.tg.tg_sig(1,k) = -0.05; % to counter lmod above...
end
if g.sys.t(k) > 10
    g.tg.tg_sig(1,k) = 0; % to counter lmod above...
end
return
