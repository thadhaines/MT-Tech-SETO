function h_sol = i_simu(k,ks,k_inc,h,bus_sim,Y_g,Y_gnc,Y_ncg,Y_nc,rec_V1,rec_V2,bo)
%Syntax: h_sol = i_simu(k,ks,k_inc,h,bus_sim,...
%                 Y_g,Y_gnc,Y_ncg,Y_nc,rec_V1,rec_V2,bo)
% 11:13 AM 18/08/97
%Purpose forms the network interface variables
% Inputs: k - the current time step
%         ks - indicates the switching times
%         k_inc - the number of time seps between switching points
%         h vector of time steps
%         bus_sim value of bus matrix at this switching time
%         Y_g - reduced Y matrix for generators
%         Y_gnc - mutual reduced Y generators-nc loads
%         Y_ncg - mutual reduced Y nc loads generators
%         Y_nc - reduced Y matrix nc loads
%         rec_V1 - voltage recovery matrix generators
%         rec_V2 - voltage recovery matrix nc loads
%         bo bus order for this switching time
% Output: h_sol - the time step at this value of ks
% Called by: s_simu
% Version 1.1
% Date:   August 1997
% Modification: add induction generator
% Version: 1.0
% Date:    March 1997
% Author:   Graham Rogers
% Copyright (c) Joe Chow All Rights Reserved

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
% % global i_dc Vdc alpha gamma dcc_pot i_dcr  i_dci Pdc
% % global r_idx i_idx ac_bus rec_ac_bus inv_ac_bus n_conv

global g

flag =1;
kdc=10*(k-1)+1;

if isempty(g.ind.n_mot)
    g.ind.n_mot = 0;
end

if isempty(g.igen.n_ig)
    g.igen.n_ig =0;
end

jay = sqrt(-1);
psi = g.sys.psi_re(:,k) + jay*g.sys.psi_im(:,k);
vmp = g.ind.vdp(:,k) + jay*g.ind.vqp(:,k);
vmpig = g.igen.vdpig(:,k) + jay*g.igen.vqpig(:,k);

if (g.ind.n_mot~=0 && g.igen.n_ig==0) 
    ntot = g.mac.n_mac + g.ind.n_mot;
    ngm = g.mac.n_mac + g.ind.n_mot;
    int_volt=[psi; vmp]; % internal voltages of generators and motors
elseif (g.ind.n_mot==0 && g.igen.n_ig~=0)
    ntot = g.mac.n_mac + g.igen.n_ig;
    ngm = g.mac.n_mac;
    int_volt=[psi; vmpig]; % internal voltages of generators and ind. generators
elseif  (g.ind.n_mot~=0 && g.igen.n_ig~=0)
    ntot = g.mac.n_mac + g.ind.n_mot + g.igen.n_ig;
    ngm = g.mac.n_mac + g.ind.n_mot;
    int_volt=[psi; vmp; vmpig]; % internal voltages of generators, motors & ind. generators
else
    int_volt = psi;
end
h_sol = h(ks);
nbus = length(bus_sim(:,1));
cur = Y_g*int_volt; % network solution currents into generators
b_v(bo(g.ncl.nload+1:nbus),1) = rec_V1*int_volt;   % bus voltage reconstruction

if g.ncl.nload~=0
    if k~=1
        kk = k-1;
    else
        kk=k;
    end
    vnc = g.sys.bus_v(g.sys.bus_int(g.ncl.load_con(:,1)),kk);% initial value
    vnc = nc_load(bus_sim,flag,Y_nc,Y_ncg,int_volt,vnc,1e-5,k,kdc);
    
    % set nc load voltages
    b_v(bo(1:g.ncl.nload),1)=vnc;
    b_v(bo(g.ncl.nload+1:nbus),1) = b_v(bo(g.ncl.nload+1:nbus),1)+rec_V2*vnc;
    cur = cur + Y_gnc*vnc;% modify generator currents for nc loads
end
% note: the dc bus voltages are the equivalent HT bus voltages
%       and not the LT bus voltages
g.sys.bus_v(g.sys.bus_int(bus_sim(:,1)),k) = b_v;
g.sys.theta(g.sys.bus_int(bus_sim(:,1)),k) = angle(b_v);
g.sys.cur_re(:,k) = real(cur(1:g.mac.n_mac));
g.sys.cur_im(:,k) = imag(cur(1:g.mac.n_mac)); % generator currents
if g.ind.n_mot~=0
    g.ind.idmot(:,k) = -real(cur(g.mac.n_mac+1:ngm));%induction motor currents
    g.ind.iqmot(:,k) = -imag(cur(g.mac.n_mac+1:ngm));%current out of network
    g.ind.s_mot(:,k) = g.sys.bus_v(g.sys.bus_int(g.ind.ind_con(:,2)),k).*(g.ind.idmot(:,k)-jay*g.ind.iqmot(:,k));
    g.ind.p_mot(:,k) = real(g.ind.s_mot(:,k));
    g.ind.q_mot(:,k) = imag(g.ind.s_mot(:,k));
end
if g.igen.n_ig~=0
    g.igen.idig(:,k) = -real(cur(ngm+1:ntot));%induction generator currents
    g.igen.iqig(:,k) = -imag(cur(ngm+1:ntot));%current out of network
    g.igen.s_igen(:,k) = g.sys.bus_v(g.sys.bus_int(g.igen.igen_con(:,2)),k)...
        .*(g.igen.idig(:,k)-jay*g.igen.iqig(:,k));
    g.igen.pig(:,k) = real(g.igen.s_igen(:,k));
    g.igen.qig(:,k) = imag(g.igen.s_igen(:,k));
end

if g.dc.n_conv ~=0
    % calculate dc voltage and current
    V0(g.dc.r_idx,1) = abs(g.sys.bus_v(g.dc.rec_ac_bus,k)).*g.dc.dcc_pot(:,7);
    V0(g.dc.i_idx,1) = abs(g.sys.bus_v(g.dc.inv_ac_bus,k)).*g.dc.dcc_pot(:,8);
    g.dc.Vdc(g.dc.r_idx,kdc) = V0(g.dc.r_idx,1).*cos(g.dc.alpha(:,kdc)) - g.dc.i_dcr(:,kdc).*g.dc.dcc_pot(:,3);
    g.dc.Vdc(g.dc.i_idx,kdc) = V0(g.dc.i_idx,1).*cos(g.dc.gamma(:,kdc)) - g.dc.i_dci(:,kdc).*g.dc.dcc_pot(:,5);
    g.dc.i_dc(g.dc.r_idx,kdc) = g.dc.i_dcr(:,kdc);
    g.dc.i_dc(g.dc.i_idx,kdc) = g.dc.i_dci(:,kdc);
end
