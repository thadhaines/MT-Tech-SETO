function i_ac = dc_cur(V,k,kdc)
%Syntax: i_ac = dc_cur(V,k,kdc)
%
%Purpose : calculates the ac current injected into the network
%          as a function of the equivalent HT voltage 
%Input:    V the complex equivalent HT voltage
%          k time step indicator, kdc is the dc time step indicator
%Output:   i_ac the pu current injection at the HT bus
%Called by nc_load
%Version 1.0
%Date:   March 1997
%Copyright Joe Chow 1991-1997 All Rights Reserved
%     %% HVDC link variables - 63
%     global  dcsp_con  dcl_con  dcc_con
%     global  r_idx  i_idx n_dcl  n_conv  ac_bus rec_ac_bus  inv_ac_bus
%     global  inv_ac_line  rec_ac_line ac_line dcli_idx
%     global  tap tapr tapi tmax tmin tstep tmaxr tmaxi tminr tmini tstepr tstepi
%     global  Vdc  i_dc P_dc i_dcinj dc_pot alpha gamma VHT dc_sig  cur_ord dcr_dsig dci_dsig
%     global  ric_idx  rpc_idx Vdc_ref dcc_pot
%     global  no_cap_idx  cap_idx  no_ind_idx  l_no_cap  l_cap
%     global  ndcr_ud ndci_ud dcrud_idx dciud_idx dcrd_sig dcid_sig
% 
%     % States
%     %line
%     global i_dcr i_dci  v_dcc
%     global di_dcr  di_dci  dv_dcc
%     global dc_dsig % added 07/13/20 -thad
%     %rectifier
%     global v_conr dv_conr
%     %inverter
%     global v_coni dv_coni
%     
%     % added to global dc
%     global xdcr_dc dxdcr_dc xdci_dc dxdci_dc angdcr angdci t_dc
%     global dcr_dc dci_dc % damping control
%     global  ldc_idx
% 
% % 
% % global  r_idx  i_idx  dcc_pot n_conv 
% % global  i_dcr  i_dci  alpha  gamma
% % 
% % 
% Vdc NOT global....
global g
jay = sqrt(-1);
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
i_ac = (P - jay*Q)./conj(V);
return
 