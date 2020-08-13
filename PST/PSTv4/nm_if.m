function nm_if(k)
% network-machine interface
%pst_var
% copy of globals for highlighting - thad - 06/11/20

% svc variables
global  svc_con n_svc svc_idx svc_pot svcll_idx
global  svc_sig
% svc user defined damping controls
global n_dcud dcud_idx svc_dsig
%states
global B_cv B_con
%dstates
global dB_cv dB_con

% tcsc variables
global  tcsc_con n_tcsc tcsvf_idx tcsct_idx
global  B_tcsc dB_tcsc
global  tcsc_sig tcsc_dsig
global  n_tcscud dtcscud_idx  %user defined damping controls

% DeltaP/omega filter variables
global  dpw_con dpw_out dpw_pot dpw_pss_idx dpw_mb_idx dpw_idx n_dpw dpw_Td_idx dpw_Tz_idx
global  sdpw1 sdpw2 sdpw3 sdpw4 sdpw5 sdpw6
global  dsdpw1 dsdpw2 dsdpw3 dsdpw4 dsdpw5 dsdpw6

%HVDC link variables
global  dcsp_con  dcl_con  dcc_con
global  r_idx  i_idx n_dcl  n_conv  ac_bus rec_ac_bus  inv_ac_bus
global  inv_ac_line  rec_ac_line ac_line dcli_idx
global  tap tapr tapi tmax tmin tstep tmaxr tmaxi tminr tmini tstepr tstepi
global  Vdc  i_dc P_dc i_dcinj dc_pot alpha gamma VHT dc_sig  cur_ord dcr_dsig dci_dsig
global  ric_idx  rpc_idx Vdc_ref dcc_pot
global  no_cap_idx  cap_idx  no_ind_idx  l_no_cap  l_cap
global  ndcr_ud ndci_ud dcrud_idx dciud_idx dcrd_sig dcid_sig

% States
%line
global i_dcr i_dci  v_dcc
global di_dcr  di_dci  dv_dcc
%rectifier
global v_conr dv_conr
%inverter
global v_coni dv_coni
% simulation control
global sw_con  scr_con

% pss design
global ibus_con  netg_con  stab_con

global g


flag=1;
mac_ind(0,k,bus_sim,flag);
mac_igen(0,k,bus_sim,flag);
mac_sub(0,k,bus_sim,flag);
mac_tra(0,k,bus_sim,flag);
mac_em(0,k,bus_sim,flag);
mdc_sig(t(k),k);
dc_cont(0,k,10*(k-1)+1,bus_sim,flag);

% Calculate current injections and bus voltages and angles
if k >= sum(k_inc(1:3))+1
    % fault cleared
    line_sim = line_pf2;
    bus_sim = bus_pf2;
    bus_int = bus_intpf2;
    Y1 = Y_gpf2;
    Y2 = Y_gncpf2;
    Y3 = Y_ncgpf2;
    Y4 = Y_ncpf2;
    Vr1 = V_rgpf2;
    Vr2 = V_rncpf2;
    bo = bopf2;
    h_sol = i_simu(k,ks,k_inc,h,bus_sim,Y1,Y2,Y3,Y4,Vr1,Vr2,bo);
elseif k >=sum(k_inc(1:2))+1
    % near bus cleared
    line_sim = line_pf1;
    bus_sim = bus_pf1;
    bus_int = bus_intpf1;
    Y1 = Y_gpf1;
    Y2 = Y_gncpf1;
    Y3 = Y_ncgpf1;
    Y4 = Y_ncpf1;
    Vr1 = V_rgpf1;
    Vr2 = V_rncpf1;
    bo = bopf1;
    h_sol = i_simu(k,ks,k_inc,h,bus_sim,Y1,Y2,Y3,Y4,Vr1,Vr2,bo);
elseif k>=k_inc(1)+1
    % fault applied
    line_sim = line_f;
    bus_sim = bus_f;
    bus_int = bus_intf;
    Y1 = Y_gf;
    Y2 = Y_gncf;
    Y3 = Y_ncgf;
    Y4 = Y_ncf;
    Vr1 = V_rgf;
    Vr2 = V_rncf;
    bo = bof;
    h_sol = i_simu(k,ks,k_inc,h,bus_sim,Y1,Y2,Y3,Y4,Vr1,Vr2,bo);
elseif k<k_inc(1)+1
    % pre fault
    line_sim = line;
    bus_sim = bus;
    bus_int = bus_intprf;
    Y1 = Y_gprf;
    Y2 = Y_gncprf;
    Y3 = Y_ncgprf;
    Y4 = Y_ncprf;
    Vr1 = V_rgprf;
    Vr2 = V_rncprf;
    bo = boprf;
    h_sol = i_simu(k,ks,k_inc,h,bus_sim,Y1,Y2,Y3,Y4,Vr1,Vr2,bo);
end
% HVDC
if ndcr_ud~=0
    % calculate the new value of bus angles rectifier user defined control
    tot_states=0;
    for jj = 1:ndcr_ud
        b_num1 = dcr_dc{jj,3};b_num2 = dcr_dc{jj,4};conv_num = dcr_dc{jj,2};
        angdcr(jj,k)=(theta(bus_int(b_num1),k)-theta(bus_int(b_num2),k));
        dcrd_sig(jj,k)=angdcr(jj,k);
        st_state = tot_states+1; dcr_states = dcr_dc{jj,7}; tot_states = tot_states+dcr_states;
        ydcrmx=dcr_dc{jj,5};ydcrmn = dcr_dc{jj,6};
        dcr_dsig(jj,k) = ...
            dcr_sud(jj,k,flag,dcr_dc{jj,1},dcrd_sig(jj,k),ydcrmx,ydcrmn,xdcr_dc(st_state:tot_states,10*(k-1)+1));
    end
end
if ndci_ud~=0
    % calculate the new value of bus angles inverter user defined control
    for jj = 1:ndci_ud
        tot_states=0;
        b_num1 = dci_dc{jj,3};b_num2 = dci_dc{jj,4};conv_num = dci_dc{jj,2};
        angdci(jj,k)=theta(bus_int(b_num1),k)-theta(bus_int(b_num2),k);
        dcid_sig(jj,k)=(angdci(jj,k)-angdci(jj,k-1))/(t(k)-t(k-1));
        st_state = tot_states+1; dci_states = dci_dc{jj,7}; tot_states = tot_states+dci_states;
        ydcimx=dci_dc{jj,5};ydcimn = dci_dc{jj,6};
        dci_dsig(jj,k) = ...
            dci_sud(jj,k,flag,dci_dc{jj,1},dcid_sig(jj,k),ydcimx,ydcimn,xdci_dc(st_state:tot_states,10*(k-1)+1));
    end
end
dc_cont(0,k,10*(k-1)+1,bus_sim,flag);
% network interface for control models
dpwf(0,k,bus_sim,flag);
pss(0,k,bus_sim,flag);

mexc_sig(k);
smpexc(0,k,flag);
exc_st3(0,k,flag);
exc_dc12(0,k,flag);

mtg_sig(t(k),k);
tg(0,k,bus_sim,flag);
tg_hydro(0,k,bus_sim,flag);
if n_dcud~=0
    % set the new line currents
    for jj=1:n_dcud
        l_num = svc_dc{jj,3};svc_num = svc_dc{jj,2};
        from_bus = bus_int(line_sim(l_num,1)); to_bus=bus_int(line_sim(l_num,2));
        svc_bn = bus_int(svc_con(svc_num,2));
        V1 = bus_v(from_bus,k);
        V2 = bus_v(to_bus,k);
        R = line_sim(l_num,3);X=line_sim(l_num,4);B=line_sim(l_num,5);tap = line_sim(l_num,6);phi = line_sim(l_num,7);
        [l_if,l_it] = line_cur(V1,V2,R,X,B,tap,phi);
        if svc_bn == from_bus; d_sig(jj,k)=abs(l_if);elseif svc_bn==to_bus;d_sig(jj,k)=abs(l_it); end
    end
end
if n_tcscud~=0
    % set the new bus voltages
    for jj=1:n_tcscud
        b_num = tcsc_dc{jj,3};tcsc_num = tcsc_dc{jj,2};
        td_sig(jj,k)=abs(bus_v(bus_int(b_num),k));
    end
end
