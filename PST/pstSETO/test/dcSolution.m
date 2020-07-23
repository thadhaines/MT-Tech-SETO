function dcSolution(k)
%DCSOLUTION Performs the dc solution for index k
% DCSOLUTION Performs the dc solution for index k
%
% Syntax: dcSolution(k)
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
%   07/23/20    11:34   Thad Haines     Version 1

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
       %% integrate dc at ten times rate (DC Stuff? 5/14/20)
        mdc_sig(k);
        if g.dc.n_conv~=0
            hdc_sol = g.k.h_sol/10;
            for kk = 1:10
                kdc=10*(k-1)+kk;
                [g.dc.xdcr_dc(:,kdc:kdc+1),g.dc.dxdcr_dc(:,kdc:kdc+1),g.dc.xdci_dc(:,kdc:kdc+1),g.dc.dxdci_dc(:,kdc:kdc+1)] = ...
                    dc_sim(k,kk,g.dc.dcr_dc,g.dc.dci_dc,g.dc.xdcr_dc(:,kdc),g.dc.xdci_dc(:,kdc),g.bus.bus_sim,hdc_sol); % dc_sim
            end
        else
            dc_cont(0,k,k,g.bus.bus_sim,2);
            dc_line(0,k,k,g.bus.bus_sim,2);
        end
end %end dcSolution