function i_ac = dc_cur(V,k,kdc)
%DC_CUR calculates the ac current injected into the network.
% DC_CUR  calculates the ac current injected into the network as a function 
% of the equivalent HT voltage.
%
% Syntax: i_ac = dc_cur(V,k,kdc)
%
%   NOTES: Called by nc_load. The Vdc variable used in this function is NOT global
% 
%   Input: 
%   V - the complex equivalent HT voltage
%   k - time step indicator
%   kdc - the dc time step indicator
%
%   Output: 
%   i_ac - the pu current injection at the HT bus
%
%   History:
%   Date        Time    Engineer        Description
%   03/xx/97    XX:XX   Joe Chow        Version 1
%   (c) Copyright Joe Chow 1991-1997 All Rights Reserved
%   07/15/20    10:42   Thad Haines     Revised format of globals and internal function documentation

% Vdc NOT global.
global g

V0(g.dc.r_idx,1) = g.dc.dcc_pot(:,7).*abs(V(g.dc.r_idx));
V0(g.dc.i_idx,1) = g.dc.dcc_pot(:,8).*abs(V(g.dc.i_idx));
dc_ang(g.dc.r_idx,1) = g.dc.alpha(:,kdc);
dc_ang(g.dc.i_idx,1) = g.dc.gamma(:,kdc);
Rc(g.dc.r_idx,1) = g.dc.dcc_pot(:,3);
Rc(g.dc.i_idx,1) = g.dc.dcc_pot(:,5);
idc(g.dc.r_idx,1) = g.dc.i_dcr(:,kdc);
idc(g.dc.i_idx,1) = g.dc.i_dci(:,kdc);
Vdc = V0.*cos(dc_ang) - idc.*Rc;
cphi = Vdc./V0;
sphi = sqrt(ones(g.dc.n_conv,1) - cphi.*cphi);
P = Vdc.*idc/g.sys.basmva;
Q = P.*sphi./cphi;
P(g.dc.i_idx) = - P(g.dc.i_idx);
i_ac = (P - 1j*Q)./conj(V);
return
 