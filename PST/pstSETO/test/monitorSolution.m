function monitorSolution(k)
%MONITORSOLUTION performs calculations related to line and area monitors
% MONITORSOLUTION performs calculations related to line and area monitors
%
% Syntax: monitorSolution(k)
%
%   NOTES:
%
%   Input:
%   k - data index to operate on
%
%   Output:
%   VOID
%
%   History:
%   Date        Time    Engineer        Description
%   07/30/20    11:13   Thad Haines     Version 1

%% Remaining 'loose' globals
% ivm variables - 5
global n_ivm mac_ivm_idx ivmmod_data ivmmod_d_sig ivmmod_e_sig

% DeltaP/omega filter variables - 21
global  dpw_con dpw_out dpw_pot dpw_pss_idx dpw_mb_idx dpw_idx n_dpw dpw_Td_idx dpw_Tz_idx
global  sdpw1 sdpw2 sdpw3 sdpw4 sdpw5 sdpw6
global  dsdpw1 dsdpw2 dsdpw3 dsdpw4 dsdpw5 dsdpw6

% pss design - 3 - Not used in Simulation? - thad 07/18/20
global ibus_con  netg_con  stab_con

%%
global g

%% ========================================================================
%% Line Monitoring and Area Calculations ==================================
%% Line Monitoring
if g.lmon.n_lmon~=0
    lmon(k)
end

%% Average Frequency Calculation
calcAveF(k,1);

%% Area Total Calcvulations
if g.area.n_area ~= 0
    calcAreaVals(k,1);
end
end % end monitorSolution