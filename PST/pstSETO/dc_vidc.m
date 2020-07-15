function dc_vidc(k,kdc)
%DC_VIDC updates Vdc and i_dc assuming ac bus voltage remains constant
% DC_VIDC  updates Vdc and i_dc assuming ac bus voltage remains constant
%
% Syntax: updates Vdc and i_dc assuming ac bus voltage remains constant
%
%   NOTES:  
% 
%   Input: 
%   k - integer time (data index)
%   kdc - integer time for dc (dc data index)
%
%   Output: 
%   VOID
%
%   History:
%   Date        Time    Engineer        Description
%   xx/xx/xx    XX:XX   xxx             Version 1.0
%   07/15/20    10:40   Thad Haines     Revised format of globals and internal function documentation

global g

V0(g.dc.r_idx,1) = abs(g.sys.bus_v(g.dc.rec_ac_bus,k)).*g.dc.dcc_pot(:,7);
V0(g.dc.i_idx,1) = abs(g.sys.bus_v(g.dc.inv_ac_bus,k)).*g.dc.dcc_pot(:,8);
g.dc.Vdc(g.dc.r_idx,kdc) = V0(g.dc.r_idx,1).*cos(g.dc.alpha(:,kdc)) - g.dc.i_dcr(:,kdc).*g.dc.dcc_pot(:,3);
g.dc.Vdc(g.dc.i_idx,kdc) = V0(g.dc.i_idx,1).*cos(g.dc.gamma(:,kdc)) - g.dc.i_dci(:,kdc).*g.dc.dcc_pot(:,5);
g.dc.i_dc(g.dc.r_idx,kdc) = g.dc.i_dcr(:,kdc);
g.dc.i_dc(g.dc.i_idx,kdc) = g.dc.i_dci(:,kdc);
