function reInitSub(i,kT)
%REINITSUB  re-initialize sub transient generator
% REINITSUB  re-initialize sub transient generator
%
% Syntax: reInitSub(i,kT)
%
%   NOTES: mostly modified code from 0 flag of mac_sub
%			Speed that is chosen to init to taken from gen 1 - should be rethought
%
%   Input:
%   i - machine number to re-init
%	kT - data index to re-int to
%
%   Output:
%   VOID
%
%   History:
%   Date        Time    Engineer        Description
%   08/28/20    14:12   Thad Haines     Version 1

global g

busnum = g.mac.mac_con(i,2);

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
% below is the non mac_pot section of the non-vectorized init from mac_sub

% extract bus information

% terminal bus voltage
g.mac.eterm(i,kT) = abs(g.bus.bus_v(conBus,kT)); % get from connected bus
% terminal bus angle in radians
g.bus.theta(busnum,kT) = angle(g.bus.bus_v(conBus,kT)); % get from connected bus

g.mac.pelect(i,kT) = 0; % system base
% electrical output power, active
g.mac.qelect(i,kT) = 0;
% electrical output power, reactive
curr = 0;  % current magnitude on genenearor base
phi = 0; % power factor angle
v = g.mac.eterm(i,kT)*exp(1j*g.bus.theta(busnum,kT)); % complex voltage in system reference frame
curr = 0;% curr*exp(1j*(g.bus.theta(busnum,1)-phi)); % complex current in system reference frame
ei = v + (g.mac.mac_con(i,5)+1j*g.mac.mac_con(i,11))*curr;
g.mac.mac_ang(i,kT) = atan2(imag(ei),real(ei)); % machine angle (delta)
g.mac.mac_spd(i,kT) = g.mac.mac_spd(1,kT); % machine speed at steady state %NOTE: figure out more generic way to do this.
rot = 1j*exp(-1j*g.mac.mac_ang(i,kT));  % system reference frame rotation
% the following should all be zero as curr == 0
curr = curr*rot;
g.mac.curdg(i,kT) = real(curr);
g.mac.curqg(i,kT) = imag(curr);% current in Park's frame
g.mac.curd(i,kT) = real(curr)/g.mac.mac_pot(i,1);% current on system base
g.mac.curq(i,kT) = imag(curr)/g.mac.mac_pot(i,1);
mcurmag = abs(curr); % current magnitude on machine base
g.mac.pmech(i,kT) = g.mac.pelect(i,kT)*g.mac.mac_pot(i,1) + g.mac.mac_con(i,5)*(mcurmag*mcurmag);
%pmech = g.mac.pelect + losses on machine base
v = v*rot;
g.mac.ed(i,kT) = real(v);
g.mac.eq(i,kT) = imag(v);
eqra = g.mac.eq(i,kT) + g.mac.mac_con(i,5)*g.mac.curqg(i,kT);
g.mac.psidpp = eqra + g.mac.mac_con(i,8)*g.mac.curdg(i,kT);
g.mac.psikd(i,kT) = eqra + g.mac.mac_con(i,4)*g.mac.curdg(i,kT);
g.mac.eqprime(i,kT) = eqra + g.mac.mac_con(i,7)*g.mac.curdg(i,kT);
edra = -g.mac.ed(i,kT)-g.mac.mac_con(i,5)*g.mac.curdg(i,kT);
g.mac.psiqpp = edra + g.mac.mac_con(i,13)*g.mac.curqg(i,kT);
g.mac.psikq(i,kT) = edra + g.mac.mac_con(i,4)*g.mac.curqg(i,kT);
g.mac.edprime(i,kT) = edra + g.mac.mac_con(i,12)*g.mac.curqg(i,kT);
% % compute saturation % handled during first init
% inv_sat = inv([0.64 0.8 1;1 1 1;1.44 1.2 1]);
% b = [0.8 1+g.mac.mac_con(i,20) 1.2*(1+g.mac.mac_con(i,21))];
% g.mac.mac_pot(i,3) = b*inv_sat(1,:)';
% g.mac.mac_pot(i,4) = b*inv_sat(2,:)';
% g.mac.mac_pot(i,5) = b*inv_sat(3,:)';
E_Isat = g.mac.mac_pot(i,3)*g.mac.eqprime(i,kT)^2 ...
    + g.mac.mac_pot(i,4)*g.mac.eqprime(i,kT) + g.mac.mac_pot(i,5);
if g.mac.eqprime(i,1)<0.8
    E_Isat=g.mac.eqprime(i,kT);
end
g.mac.vex(i,kT) = E_Isat + g.mac.mac_pot(i,6)*(g.mac.eqprime(i,kT)-...
    g.mac.psikd(i,kT))+g.mac.mac_pot(i,7)*g.mac.curdg(i,kT);
g.mac.fldcur(i,kT) = g.mac.vex(i,kT);
g.mac.psi_re(i,kT) = sin(g.mac.mac_ang(i,kT)).*(-g.mac.psiqpp) + ...
    cos(g.mac.mac_ang(i,kT)).*g.mac.psidpp; % real part of psi
g.mac.psi_im(i,kT) = -cos(g.mac.mac_ang(i,kT)).*(-g.mac.psiqpp) + ...
    sin(g.mac.mac_ang(i,kT)).*g.mac.psidpp; % imag part of psi
end