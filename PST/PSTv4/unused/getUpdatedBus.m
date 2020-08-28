function busNew = getUpdatedBus(k)
%GETUPDATEDBUS  returns bus matrix with collected values from index k
% GETUPDATEDBUS  returns bus matrix with collected values from index k
%
% Syntax: busNew = getUpdatedBus(k)
%
%   NOTES: maybe not super useful...
%
%   Input:
%   k - data index to 
%
%   Output:
%   busNew - bus array with updated values from index k
%
%   History:
%   Date        Time    Engineer        Description
%   08/27/20    09:37   Thad Haines     Version 1

%%
global g

busNew = g.bus.bus_sim;

% Assert that the following values do not change over the simulation:
% col1 bus number
% col6 p_load(pu)
% col7 q_load(pu)
% col8 G shunt(pu)
% col9 B shunt(pu)
% col10 bus_type
% col11 q_gen_max(pu)
% col12 q_gen_min(pu)
% col13 v_rated (kV)
% col14 v_max  pu
% col15 v_min  pu

% Collect values that do change over simulation
% bus voltage
nfBusN = size(g.bus.bus_v,1)-1; % size of non-fault bus voltages
sLoc = find(g.bus.bus(:,10) == 1); % slack gen loc
busNew(:,2) = abs(g.bus.bus_v(1:nfBusN,k));

% bus angle in degree
busNew(:,3) = rad2deg(angle(g.bus.bus_v(1:nfBusN,k)));
busNew(:,3) = busNew(:,3)- busNew(sLoc,3); % adjust to be in reference to slack

% p gen from machines
macBusIndx = find(busNew(:,10) == 1 | busNew(:,10) == 2);
busNew(macBusIndx, 4) = g.mac.pelect(:,k);
% q gen from machines
busNew(macBusIndx, 5) = g.mac.qelect(:,k);

end
