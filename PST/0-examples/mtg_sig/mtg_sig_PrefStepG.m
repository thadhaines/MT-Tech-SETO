function mtg_sig(k)
% MTG_SIG Defines modulation signal for turbine power reference
% Syntax: mtg_sig(k)
%
%   History:
%   Date        Time    Engineer        Description
%   07/xx/98    12:37   Graham Rogers   Version 1
%   06/05/20    16:21   Thad Haines     V2 - using global g, no t passed
%   in, no dummy varibale f passed out. Default behavior commented out.
%

global g

if g.sys.t(k) > 0.5
    g.tg.tg_sig(1,k) = 0.05; 
end
if g.sys.t(k) > 5
    g.tg.tg_sig(1,k) = -0.05; 
end
if g.sys.t(k) > 10
    g.tg.tg_sig(1,k) = 0;
end
return
