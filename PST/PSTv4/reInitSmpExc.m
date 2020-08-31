function reInitSmpExc(exNum, kT)
%REINITSMPEXC  re-initialize simple exciter
% REINITSMPEXC  re-initialize simple exciter
%
% Syntax: reInitSmpExc(exNum, kT)
%
%   NOTES: mostly modified code from 0 flag of smpexc
%
%   Input:
%   exNum - exciter number to re-init
%	kT - data index to re-int to
%
%   Output:
%   VOID
%
%   History:
%   Date        Time    Engineer        Description
%   08/31/20    09:25   Thad Haines     Version 1

global g


n = g.mac.mac_int(g.exc.exc_con(exNum,2)); % machine number


busnum = g.mac.mac_con(n,2);
% find bus generator is re-connecting to
% machine is connected as from
fBus = find(g.line.line(:,1) == busnum);
tBus = find(g.line.line(:,2) == busnum);
if ~isempty(fBus)
    conBus = g.line.line(fBus,2);
% machine is connected as to
elseif ~isempty(tBus)
    conBus = g.line.line(fBus,1);
else
    error('bus connecting to generator not found.')
end


g.exc.Efd(exNum,kT) = g.mac.vex(n,kT); % machine field foltage
g.exc.V_A(exNum,kT) = g.exc.Efd(exNum,kT)/g.exc.exc_con(exNum,4); % laglead
g.exc.V_As(exNum,kT) = g.exc.V_A(exNum,kT); % leadlag state variable
g.exc.V_TR(exNum,kT) =   g.mac.eterm(n,kT); % input filter state % machine terminal voltage


err = g.exc.V_A(exNum,kT); % summing junction error

if g.exc.exc_con(exNum,6)~=0
    g.exc.exc_pot(exNum,5) = g.exc.exc_con(exNum,7)/g.exc.exc_con(exNum,6);
else
    g.exc.exc_pot(exNum,5)=1;
end

g.exc.exc_pot(exNum,3) = g.mac.eterm(n,kT)+err; % reference 

end