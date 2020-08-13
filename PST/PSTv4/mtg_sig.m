function mtg_sig(k)
% MTG_SIG Defines modulation signal for turbine power reference
% Syntax: f = mtg_sig(k)
%
%   NOTE: Can be totally blank IF tg_sig not manipulated in any way (i.e. AGC)
%
%   History:
%   Date        Time    Engineer        Description
%   07/xx/98    12:37   Graham Rogers   Version 1
%   06/05/20    16:21   Thad Haines     V2 - using global g, no t passed
%   in, no dummy varibale f passed out. Default behavior commented out.
%

global g
g.tg.tg_sig(:,k) = zeros(size(g.tg.tg_sig,1),1);

end% end function