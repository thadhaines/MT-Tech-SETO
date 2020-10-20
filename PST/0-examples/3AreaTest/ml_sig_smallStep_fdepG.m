function ml_sig(k)
% Syntax: ml_sig(k)
%16:45 PM 29/10/20
% defines modulation signal for lmod control

global g

% Add frequency dependant load with 1:1 sensitivity
if g.sys.t(k) == 0
    disp('Adding frequency dependent load characteristics')
end
g.lmod.lmod_sig(1,k) = (1-g.mac.mac_spd(1,k)) * g.bus.bus( g.area.area(1).loadBusNdx(1), 6);
g.lmod.lmod_sig(2,k) = (1-g.mac.mac_spd(2,k)) * g.bus.bus( g.area.area(2).loadBusNdx(1), 6);
g.lmod.lmod_sig(3,k) = (1-g.mac.mac_spd(3,k)) * g.bus.bus( g.area.area(3).loadBusNdx(1), 6);

if g.sys.t(k) > 2
    % load step
    g.lmod.lmod_sig(1,k) = g.lmod.lmod_sig(1,k) + 0.02;
end

return