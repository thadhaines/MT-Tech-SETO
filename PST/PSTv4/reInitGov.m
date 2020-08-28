function reInitGov(tgNum, kT)
%REINITGOV  re-initialize governor
% REINITGOV  re-initialize governor
%
% Syntax: reInitGov(tgNum, kT)
%
%   NOTES: mostly modified code from 0 flag of tg
%
%   Input:
%   tgNum - governor number to re-init
%	kT - data index to re-int to
%
%   Output:
%   VOID
%
%   History:
%   Date        Time    Engineer        Description
%   08/28/20    14:12   Thad Haines     Version 1

global g

macNum = g.tg.tg_con(tgNum,2); % get connected machine number 

% Check for pmech being inside generator limits
if g.mac.pmech(macNum,kT) > g.tg.tg_con(tgNum,5)
    error('TG init: pmech > upper limit, check machine base')
end
if g.mac.pmech(macNum,kT) < 0
    error('TG init: pmech < 0, check data')
end

% re-Initialize states
g.tg.tg1(tgNum,kT) = g.mac.pmech(macNum,kT);
%
g.tg.tg_pot(tgNum,1) = g.tg.tg_con(tgNum,8)/g.tg.tg_con(tgNum,7);
a1 = 1 - g.tg.tg_pot(tgNum,1);
g.tg.tg_pot(tgNum,2) = a1;
g.tg.tg2(tgNum,kT) = a1*g.mac.pmech(macNum,kT);
%
g.tg.tg_pot(tgNum,3) = g.tg.tg_con(tgNum,9)/g.tg.tg_con(tgNum,10);
a2 = 1 - g.tg.tg_pot(tgNum,3);
g.tg.tg_pot(tgNum,4) = a2;
g.tg.tg3(tgNum,kT) = a2*g.mac.pmech(macNum,kT);
%
g.tg.tg_pot(tgNum,5) = g.mac.pmech(macNum,kT); % Pref
%
g.tg.tg_sig(tgNum,kT) = 0;
% 
% % set derivatives to zero
% g.tg.dtg1(tgNum, kT) = 0;
% g.tg.dtg2(tgNum, kT) = 0;
% g.tg.dtg3(tgNum, kT) = 0;
% g.tg.dtg4(tgNum, kT) = 0;
% g.tg.dtg5(tgNum, kT) = 0;

end