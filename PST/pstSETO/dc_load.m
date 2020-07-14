function [Yrr,Yri,Yir,Yii] = dc_load(V,k,kdc)
%Syntax: [Yrr,Yri,Yir,Yii] = dc_load(V,k,kdc)
%Purpose: To calculate the non-linear Jacobian elements
%         associated with a line commutated HVDC link
% Inputs:
%         V - Equivalent HT terminal voltage
%         k - step indicator; kdc dc time step indicator
% Outputs:
%         Yrr - dir/dvr
%         Yri - dir/dvi
%         Yir - dii/dvr
%         Yii - dii/dvi
% Called by: nc_load
%Version: 1.0
%Date:    March 1997
%Author:  Graham Rogers
% Copyright (c) Joe Chow 1991-1997 All Rights Reserved

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
% % global  i_dci  i_dcr  dcc_pot  alpha  gamma  r_idx  i_idx
% % global  n_conv n_dcl

%Vdc NOT Global
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
iac = (P - jay*Q)./conj(V);
ir = real(iac);
ii = imag(iac);
dV0dVr(g.dc.r_idx,1) = g.dc.dcc_pot(:,7).*real(V(g.dc.r_idx))./abs(V(g.dc.r_idx));
dV0dVr(g.dc.i_idx,1) = g.dc.dcc_pot(:,8).*real(V(g.dc.i_idx))./abs(V(g.dc.i_idx));
dV0dVi(g.dc.r_idx,1) = g.dc.dcc_pot(:,7).*imag(V(g.dc.r_idx))./abs(V(g.dc.r_idx));
dV0dVi(g.dc.i_idx,1) = g.dc.dcc_pot(:,8).*imag(V(g.dc.i_idx))./abs(V(g.dc.i_idx));
dPdVr = idc.*cos(dc_ang).*dV0dVr/g.sys.basmva;
dPdVi = idc.*cos(dc_ang).*dV0dVi/g.sys.basmva;
Kq = idc.*(ones(g.dc.n_conv,1)-cos(dc_ang).*cphi)./sphi/g.sys.basmva;
dQdVr = Kq.*dV0dVr;
dQdVi = Kq.*dV0dVi;
Vr = real(V);
Vi = imag(V);
Vmag2 = Vr.*Vr + Vi.*Vi;
Vmag4 = Vmag2.*Vmag2;
Yrr = (P + Vi.*dQdVr + Vr.*dPdVr)./Vmag2;
Yrr = Yrr - 2*(P.*Vr + Q.*Vi).*Vr./Vmag4;
Yri = (Q + Vr.*dPdVi + Vi.*dQdVi)./Vmag2;
Yri = Yri - 2*(P.*Vr + Q.*Vi).*Vi./Vmag4;
Yir = -(Q - Vi.*dPdVr + Vr.*dQdVr)./Vmag2;
Yir = Yir + 2*(Q.*Vr - P.*Vi).*Vr./Vmag4;
Yii = (P + Vi.*dPdVi - Vr.*dQdVi)./Vmag2;
Yiii = Yii + 2*(Q.*Vr - P.*Vi).*Vi./Vmag4;
