function initStep(k)
%INITSTEP Performs intialization performed at the beginning of a solution step
% INITSTEP Performs intialization performed at the beginning of a solution step
% functions lifed from s_simu_Batch
%
% Syntax: initStep(k)
%
%   NOTES:
%
%   Input:
%   k - data index
%
%   Output:
%   VOID
%
%   History:
%   Date        Time    Engineer        Description
%   07/27/20    10:27   Thad Haines     Version 1

global g

g.sys.mach_ref(k) = 0;
g.mac.pmech(:,k+1) = g.mac.pmech(:,k);
g.igen.tmig(:,k+1) = g.igen.tmig(:,k);

if g.dc.n_conv~=0
    g.dc.cur_ord(:,k+1) = g.dc.cur_ord(:,k);
end

% Trip gen - Copied from v2.3 06/01/20 - thad
[f,g.mac.mac_trip_states] = mac_trip_logic(g.mac.mac_trip_flags,g.mac.mac_trip_states,g.sys.t,k);
g.mac.mac_trip_flags = g.mac.mac_trip_flags | f;

end % end initStep
