function networkSolution(k)
%NETWORKSOLUTION Performs the network solution for index k
% NETWORKSOLUTION Performs the network solution for index k
%
% Syntax: networkSolution(k)
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
%   07/23/20    11:28   Thad Haines     Version 1
%   09/03/20    12:00   Thad Haines     Version 1.1 - changed g.int to g.y

%% Remaining 'loose' globals
% DeltaP/omega filter variables - 21
global  dpw_con dpw_out dpw_pot dpw_pss_idx dpw_mb_idx dpw_idx n_dpw dpw_Td_idx dpw_Tz_idx
global  sdpw1 sdpw2 sdpw3 sdpw4 sdpw5 sdpw6
global  dsdpw1 dsdpw2 dsdpw3 dsdpw4 dsdpw5 dsdpw6

% pss design - 3 - Not used in Simulation? - thad 07/18/20
global ibus_con  netg_con  stab_con

%%
global g

%% Flag = 1
flag = 1;
%timestep = int2str(k); % not used? 06/09/20
g.sys.mach_ref(k) = 0;
% network-machine interface
mac_ind(0,k,g.bus.bus_sim,flag);
mac_igen(0,k,g.bus.bus_sim,flag);
mac_sub(0,k,g.bus.bus_sim,flag);
mac_tra(0,k,g.bus.bus_sim,flag);
mac_em(0,k,g.bus.bus_sim,flag);
mac_ivm(0,k,g.bus.bus_sim,flag);

mdc_sig(k); % dc controls mod signals
dc_cont(0,k,10*(k-1)+1,g.bus.bus_sim,flag); % Models the action of HVDC link pole controllers

%% Swtich System based on sw_con trip settings... ======================
% ======================================================================

% Calculate current injections and bus voltages and angles
if k >= sum(g.k.k_inc(1:3))+1
    %% fault cleared - post fault 2
    g.line.line_sim = g.line.line_pf2;
    g.bus.bus_sim = g.bus.bus_pf2;
    g.bus.bus_int = g.bus.bus_intpf2;
    
    Y1 = g.y.Y_gpf2;
    Y2 = g.y.Y_gncpf2;
    Y3 = g.y.Y_ncgpf2;
    Y4 = g.y.Y_ncpf2;
    Vr1 = g.y.V_rgpf2;
    Vr2 = g.y.V_rncpf2;
    bo = g.y.bopf2;
    
    % i_simu forms the network interface variables
    %h_sol = i_simu(k,ks,k_inc,h,bus_sim,Y1,Y2,Y3,Y4,Vr1,Vr2,bo);
    % duplicate call?
    % h_sol calculated after this 'if' block...
    
elseif k >=sum(g.k.k_inc(1:2))+1
    %% near bus cleared - post fault 1
    g.line.line_sim = g.line.line_pf1;
    
    g.bus.bus_sim = g.bus.bus_pf1;
    g.bus.bus_int = g.bus.bus_intpf1;
    
    Y1 = g.y.Y_gpf1;
    Y2 = g.y.Y_gncpf1;
    Y3 = g.y.Y_ncgpf1;
    Y4 = g.y.Y_ncpf1;
    Vr1 = g.y.V_rgpf1;
    Vr2 = g.y.V_rncpf1;
    bo = g.y.bopf1;
    
    %h_sol = i_simu(k,ks,k_inc,h,bus_sim,Y1,Y2,Y3,Y4,Vr1,Vr2,bo);
    
elseif k>=g.k.k_inc(1)+1
    %% fault applied - fault
    g.line.line_sim = g.line.line_f;
    g.bus.bus_sim = g.bus.bus_f;
    
    g.bus.bus_int = g.bus.bus_intf;
    
    Y1 = g.y.Y_gf;
    Y2 = g.y.Y_gncf;
    Y3 = g.y.Y_ncgf;
    Y4 = g.y.Y_ncf;
    Vr1 = g.y.V_rgf;
    Vr2 = g.y.V_rncf;
    bo = g.y.bof;
    
    %h_sol = i_simu(k,ks,k_inc,h,bus_sim,Y1,Y2,Y3,Y4,Vr1,Vr2,bo);
    
elseif k<g.k.k_inc(1)+1
    %% pre fault
    g.line.line_sim = g.line.line;
    g.bus.bus_sim = g.bus.bus;
    
    g.bus.bus_int = g.bus.bus_intprf;
    
    Y1 = g.y.Y_gprf;
    Y2 = g.y.Y_gncprf;
    Y3 = g.y.Y_ncgprf;
    Y4 = g.y.Y_ncprf;
    Vr1 = g.y.V_rgprf;
    Vr2 = g.y.V_rncprf;
    bo = g.y.boprf;
    
    %h_sol = i_simu(k,ks,k_inc,h,bus_sim,Y1,Y2,Y3,Y4,Vr1,Vr2,bo);
end

%% apply gen trip - added from v2.3 - 06/01/20 - thad
if sum(g.mac.mac_trip_flags)>0.5
    genBuses = g.mac.mac_con(g.mac.mac_trip_flags==1,2);
    for kB=1:length(genBuses)
        nL = find(genBuses(kB)==g.line.line_sim(:,1) | genBuses(kB)==g.line.line_sim(:,2));
        if isempty(nL)
            error('*!* Line connecting generator to trip not found'); 
        end
        g.line.line_sim(nL,4) = 1e7; %make reactance infinity
    end
    [Y1,Y2,Y3,Y4,Vr1,Vr2,bo] = red_ybus(g.bus.bus_sim,g.line.line_sim);
    clear nL kB genBuses
end

%% Solve Network solution
g.k.h_sol = i_simu(k,g.k.ks,g.k.k_inc,g.k.h,g.bus.bus_sim,Y1,Y2,Y3,Y4,Vr1,Vr2,bo);

% Executed in step 2, added here for consistancy... -thad 07/23/20
g.mac.vex(:,k+1) = g.mac.vex(:,k);
g.dc.cur_ord(:,k+1) = g.dc.cur_ord(:,k);

%% HVDC
if g.dc.ndcr_ud~=0
    % calculate the new value of bus angles rectifier user defined control
    tot_states=0;
    for jj = 1:g.dc.ndcr_ud
        b_num1 = g.dc.dcr_dc{jj,3};
        b_num2 = g.dc.dcr_dc{jj,4};
        conv_num = g.dc.dcr_dc{jj,2};
        g.dc.angdcr(jj,k) = (g.bus.theta(g.bus.bus_int(b_num1),k)-g.bus.theta(g.bus.bus_int(b_num2),k));
        g.dc.dcrd_sig(jj,k)=g.dc.angdcr(jj,k);
        st_state = tot_states+1;
        dcr_states = g.dc.dcr_dc{jj,7};
        tot_states = tot_states+dcr_states;
        ydcrmx= g.dc.dcr_dc{jj,5};
        ydcrmn = g.dc.dcr_dc{jj,6};
        g.dc.dcr_dsig(jj,k) = ...
            dcr_sud(jj,k,flag,g.dc.dcr_dc{jj,1},g.dc.dcrd_sig(jj,k),ydcrmx,ydcrmn,g.dc.xdcr_dc(st_state:tot_states,10*(k-1)+1));
    end
end
if g.dc.ndci_ud~=0
    % calculate the new value of bus angles inverter user defined control
    for jj = 1:g.dc.ndci_ud
        tot_states=0;
        b_num1 = g.dc.dci_dc{jj,3};
        b_num2 = g.dc.dci_dc{jj,4};
        conv_num = g.dc.dci_dc{jj,2};
        g.dc.angdci(jj,k)=g.bus.theta(g.bus.bus_int(b_num1),k)-g.bus.theta(g.bus.bus_int(b_num2),k);
        g.dc.dcid_sig(jj,k)=(g.dc.angdci(jj,k)-g.dc.angdci(jj,k-1))/(t(k)-t(k-1));
        st_state = tot_states+1;
        dci_states = g.dc.dci_dc{jj,7};
        tot_states = tot_states+dci_states;
        ydcimx= g.dc.dci_dc{jj,5};
        ydcimn = g.dc.dci_dc{jj,6};
        g.dc.dci_dsig(jj,k) = ...
            dci_sud(jj,k,flag,g.dc.dci_dc{jj,1},g.dc.dcid_sig(jj,k),ydcimx,ydcimn,g.dc.xdci_dc(st_state:tot_states,10*(k-1)+1));
    end
end

dc_cont(0,k,10*(k-1)+1,g.bus.bus_sim,flag);

%% network interface for control models
dpwf(0,k,flag);

mexc_sig(k);  % executed twice? - thad 07/18/20
smpexc(0,k,flag);
smppi(0,k,flag);
exc_st3(0,k,flag);
exc_dc12(0,k,flag);

mtg_sig(k); % executed twice? - thad 07/18/20
tg(0,k,flag);
tg_hydro(0,k,g.bus.bus_sim,flag);

if g.svc.n_dcud~=0
    %% set the new line currents
    % SVC damping control...
    % Static Var Compensator
    for jj=1:g.svc.n_dcud
        l_num = g.svc.svc_dc{jj,3};
        svc_num = g.svc.svc_dc{jj,2};
        from_bus = g.bus.bus_int(line_sim(l_num,1));
        to_bus = g.bus.bus_int(line_sim(l_num,2));
        svc_bn = g.bus.bus_int(g.svc.svc_con(svc_num,2));
        V1 = g.bus.bus_v(from_bus,k);
        V2 = g.bus.bus_v(to_bus,k);
        R = line_sim(l_num,3);
        X = line_sim(l_num,4);
        B = line_sim(l_num,5);
        g.dc.tap = line_sim(l_num,6);
        phi = line_sim(l_num,7);
        [l_if,l_it] = line_cur(V1,V2,R,X,B,g.dc.tap,phi);
        
        if svc_bn == from_bus
            d_sig(jj,k) = abs(l_if);
        elseif svc_bn == to_bus
            d_sig(jj,k) = abs(l_it);
        end
    end
end

if g.tcsc.n_tcscud~=0
    % set the new bus voltages
    % tcsc damping controls
    % Thyristor Controlled Series compensator
    for jj=1:g.tcsc.n_tcscud
        b_num = g.tcsc.tcsc_dc{jj,3};
        tcsc_num = g.tcsc.tcsc_dc{jj,2};
        g.tcsc.td_sig(jj,k)=abs(g.bus.bus_v(g.bus.bus_int(b_num),k));
    end
end


end% end networkSolution