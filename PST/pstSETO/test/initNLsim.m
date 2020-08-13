function initNLsim()
%INITNLSIM performs intialization for non-linear simulation
% INITNLSIM performs intialization for non-linear simulation
%
% Syntax: initNLsim()
%
%   NOTES:
%
%   Input:
%   VOID
%
%   Output:
%   VOID
%
%   History:
%   Date        Time    Engineer        Description
%   07/30/20    10:54   Thad Haines     Version 1
%   08/05/20    13:15   Thad Haines     Version 1.0.1 - added pwrmod signals to global
%   08/11/20    11:25   Thad Haines     Version 1.0.2 - added ivmmod

%% Remaining 'loose' globals
% DeltaP/omega filter variables - 21
global  dpw_con dpw_out dpw_pot dpw_pss_idx dpw_mb_idx dpw_idx n_dpw dpw_Td_idx dpw_Tz_idx
global  sdpw1 sdpw2 sdpw3 sdpw4 sdpw5 sdpw6
global  dsdpw1 dsdpw2 dsdpw3 dsdpw4 dsdpw5 dsdpw6

% pss design - 3 - Not used in Simulation? - thad 07/18/20
global ibus_con  netg_con  stab_con

%%
global g

%% step 1: construct reduced Y matrices
warning('*** Initialize Y matrix (matracies?) and Dynamic Models')
disp('constructing reduced y matrices')
disp('initializing motor,induction generator, svc and dc control models')

g.bus.bus = mac_ind(0,1,g.bus.bus,0);% initialize induction motor
g.bus.bus = mac_igen(0,1,g.bus.bus,0); % initialize induction generator
g.bus.bus = svc(0,1,g.bus.bus,0);%initialize svc
dc_cont(0,1,1,g.bus.bus,0);% initialize dc controls
% this has to be done before red_ybus is used since the motor and svc
% initialization alters the bus matrix and dc parameters are required

y_switch % calculates the reduced y matrices for the different switching conditions

disp('initializing other models...')

%% step 2: initialization
n_bus = length(g.bus.bus(:,1));
g.bus.theta(1:n_bus,1) = g.bus.bus(:,3)*pi/180;
g.bus.bus_v(1:n_bus,1) = g.bus.bus(:,2).*exp(1j*g.bus.theta(1:n_bus,1)); %complex bus voltage

% clear temp variables
clear z zdc zdcl ze zig zm t n_bus

if g.svc.n_dcud ~=0 % Seems like this should be put in a seperate script - thad 06/08/20
    %% calculate the initial magnitude of line current for svc damping controls
    for j=1:g.svc.n_dcud
        l_num = g.svc.svc_dc{j,3};
        svc_num = g.svc.svc_dc{j,2};
        from_bus = g.bus.bus_int(g.line.line(l_num,1));
        to_bus = g.bus.bus_int(g.line.line(l_num,2));
        svc_bn = g.bus.bus_int(g.svc.svc_con(svc_num,2));
        
        if svc_bn~= from_bus&& svc_bn  ~= to_bus
            error('the svc is not at the end of the specified line');
        end
        
        V1 = g.bus.bus_v(from_bus,1);
        V2 = g.bus.bus_v(to_bus,1);
        R = g.line.line(l_num,3);
        X = g.line.line(l_num,4);
        B = g.line.line(l_num,5);
        g.dc.tap = line(l_num,6);
        phi = g.line.line(l_num,7);
        [l_if,l_it] = line_cur(V1,V2,R,X,B,g.dc.tap,phi);
        l_if0(j)=l_if;
        l_it0(j)=l_it;
        
        if svc_bn == from_bus
            d_sig(j,1) = abs(l_if);
        elseif svc_bn == to_bus
            d_sig(j,1) = abs(l_it);
        end
    end
    clear j
end

if g.tcsc.n_tcscud ~=0 % Seems like this should be put in a seperate script - thad 06/08/20
    %% calculate the initial magnitude of bus voltage magnitude for tcsc damping controls
    for j=1:g.tcsc.n_tcscud
        b_num = g.tcsc.tcsc_dc{j,3};
        tcsc_num = g.tcsc.tcsc_dc{j,2};
        g.tcsc.td_sig(j,1) =abs(g.bus.bus_v(g.bus.bus_int(b_num),1));
    end
    clear j
end

if g.dc.n_conv~=0 % Seems like this should be put in a seperate script - thad 06/08/20
    %% change dc buses from LT to HT
    Pr = g.bus.bus(g.dc.rec_ac_bus,6);
    Pi = g.bus.bus(g.dc.inv_ac_bus,6);
    Qr = g.bus.bus(g.dc.rec_ac_bus,7);
    Qi = g.bus.bus(g.dc.inv_ac_bus,7);
    VLT= g.bus.bus_v(g.dc.ac_bus,1);
    i_acr = (Pr-1j*Qr)./conj(VLT(g.dc.r_idx));
    i_aci = (Pi - 1j*Qi)./conj(VLT(g.dc.i_idx));
    IHT(g.dc.r_idx,1)=i_acr;
    IHT(g.dc.i_idx,1)=i_aci;
    g.dc.VHT(g.dc.r_idx,1) = (VLT(g.dc.r_idx) + 1j*g.dc.dcc_pot(:,2).*i_acr);
    g.dc.VHT(g.dc.i_idx,1) = (VLT(g.dc.i_idx) + 1j*g.dc.dcc_pot(:,4).*i_aci);
    g.bus.bus_v(g.dc.ac_bus,1) = g.dc.VHT;
    g.bus.theta(g.dc.ac_bus,1) = angle(g.bus.bus_v(g.dc.ac_bus,1));
    % modify the bus matrix to the HT buses
    g.bus.bus(g.dc.ac_bus,2) = abs(g.bus.bus_v(g.dc.ac_bus,1));
    g.bus.bus(g.dc.ac_bus,3) = g.bus.theta(g.dc.ac_bus,1)*180/pi;
    SHT = g.dc.VHT.*conj(IHT);
    g.bus.bus(g.dc.ac_bus,6) = real(SHT);
    g.bus.bus(g.dc.ac_bus,7) = imag(SHT);
    
    if g.dc.ndcr_ud~=0 % Seems like this should be put in a seperate script - thad 06/08/20
        % calculate the initial value of bus angles rectifier user defined control
        for j = 1:g.dc.ndcr_ud
            b_num1 = g.dc.dcr_dc{j,3};
            b_num2 = g.dc.dcr_dc{j,4};
            conv_num = g.dc.dcr_dc{j,2};
            g.dc.angdcr(j,:) = g.bus.theta(g.bus.bus_int(b_num1),1)-g.bus.theta(g.bus.bus_int(b_num2),1);
            g.dc.dcrd_sig(j,:)=g.dc.angdcr(j,:);
        end
        clear j
    end
    if g.dc.ndci_ud~=0 % Seems like this should be put in a seperate script - thad 06/08/20
        % calculate the initial value of bus angles inverter user defined control
        for j = 1:g.dc.ndci_ud
            b_num1 = g.dc.dci_dc{j,3};
            b_num2 = g.dc.dci_dc{j,4};
            conv_num = g.dc.dci_dc{j,2};
            g.dc.angdci(j,:) = g.bus.theta(g.bus.bus_int(b_num1),1)-g.bus.theta(g.bus.bus_int(b_num2),1);
            g.dc.dcid_sig(j,:) = g.dc.angdci(j,:);
        end
        clear j
    end
end

% Flag = 0 == Initialization
warning('*** Dynamic model initialization via functions/scripts:')
flag = 0;
g.bus.bus_int = g.bus.bus_intprf;% pre-fault system

disp('generators')
mac_sub(0,1,g.bus.bus,flag); % first
mac_tra(0,1,g.bus.bus,flag);
mac_em(0,1,g.bus.bus,flag);
mac_ivm(0,1,g.bus.bus,flag); % ivm - thad 06/01/20

disp('generator controls')
dpwf(0,1,flag);
pss(0,1,flag);

% exciters
smpexc(0,1,flag);
smppi(0,1,flag);
exc_st3(0,1,flag);
exc_dc12(0,1,flag);

tg(0,1,flag); % modified 06/05/20 to global g
tg_hydro(0,1,g.bus.bus,flag);

%% initialize ivm modulation control - added from v2.3 06/01/20 - thad
% Seems like this should be put in a seperate script - thad 06/08/20
if g.ivm.n_ivm~=0
    disp('ivm modulation')
    [~,~,~,~,Dini,Eini] = ivmmod_dyn([],[],g.bus.bus,g.sys.t,1,flag);
    
    if (~iscell(Dini) || ~iscell(Eini))
        error('Error in ivmmod_dyn, initial states must be cells');
    end
    if (any(size(Dini)-[g.ivm.n_ivm 1]) || any(size(Eini)-[g.ivm.n_ivm 1]))
        error('Dimension error in ivmmod_dyn');
    end
    
    g.ivm.ivmmod_d_sigst = cell(g.ivm.n_ivm,1);
    g.ivm.ivmmod_e_sigst = g.ivm.ivmmod_d_sigst;
    g.ivm.divmmod_d_sigst = g.ivm.ivmmod_d_sigst;
    g.ivm.divmmod_e_sigst = g.ivm.ivmmod_d_sigst;
    
    k = size(g.sys.t,2); % length of logged values
    
    for index=1:g.ivm.n_ivm
        if ((size(Dini{index},2)~=1) || (size(Eini{index},2)~=1))
            error('Dimension error in ivmmod_dyn');
        end
        g.ivm.divmmod_d_sigst{index} = zeros(size(Dini{index},1),k);
        g.ivm.ivmmod_d_sigst{index} = Dini{index}*ones(1,k);
        g.ivm.divmmod_e_sigst{index} = zeros(size(Eini{index},1),k);
        g.ivm.ivmmod_e_sigst{index} = Eini{index}*ones(1,k);
    end
    clear index Dini Eini
end

%% initialize svc damping controls
% Seems like this should be put in a seperate script - thad 06/08/20
if g.svc.n_dcud~=0
    disp('svc damping controls')
    tot_states=0;
    for i = 1:g.svc.n_dcud
        ysvcmx = g.svc.svc_dc{i,4};
        ysvcmn = g.svc.svc_dc{i,5};
        svc_num = g.svc.svc_dc{i,2};
        st_state = tot_states+1;
        svc_states = g.svc.svc_dc{i,6};
        tot_states = tot_states+svc_states;
        [g.svc.svc_dsig(svc_num,1),g.svc.xsvc_dc(st_state:tot_states,1),g.svc.dxsvc_dc(st_state:tot_states,1)] =...
            svc_sud(i,1,flag,g.svc.svc_dc{i,1},d_sig(i,1),ysvcmx,ysvcmn);
    end
    clear i
end

%% initialize tcsc damping controls
% Seems like this should be put in a seperate script - thad 06/08/20
if g.tcsc.n_tcscud~=0
    disp('tcsc damping controls')
    tot_states=0;
    for i = 1:g.tcsc.n_tcscud
        ytcscmx = g.tcsc.tcsc_dc{i,4};
        ytcscmn = g.tcsc.tcsc_dc{i,5};
        tcsc_num = g.tcsc.tcsc_dc{i,2};
        st_state = tot_states+1;
        tcsc_states = g.tcsc.tcsc_dc{i,6};
        tot_states = tot_states+tcsc_states;
        [g.tcsc.tcsc_dsig(tcsc_num,1),g.tcsc.xtcsc_dc(st_state:tot_states,1),g.tcsc.dxtcsc_dc(st_state:tot_states,1)] =...
            tcsc_sud(i,1,flag,g.tcsc.tcsc_dc{i,1},g.tcsc.td_sig(i,1),ytcscmx,ytcscmn);
    end
    clear i
end

%% initialize rectifier damping controls
% Seems like this should be put in a seperate script - thad 06/08/20
if g.dc.ndcr_ud~=0
    disp('rectifier damping controls')
    tot_states=0;
    for i = 1:g.dc.ndcr_ud
        ydcrmx = g.dc.dcr_dc{i,5};
        ydcrmn = g.dc.dcr_dc{i,6};
        rec_num = g.dc.dcr_dc{i,2};
        st_state = tot_states+1;
        dcr_states = g.dc.dcr_dc{i,7};
        tot_states = tot_states+dcr_states;
        [g.dc.dcr_dsig(rec_num,1),g.dc.xdcr_dc(st_state:tot_states,1),g.dc.dxdcr_dc(st_state:tot_states,1)] = ...
            dcr_sud(i,1,flag,g.dc.dcr_dc{i,1},g.dc.dcrd_sig(i,1),ydcrmx,ydcrmn);
    end
    clear i
end

%% initialize inverter damping controls
% Seems like this should be put in a seperate script - thad 06/08/20
if g.dc.ndci_ud~=0
    disp('inverter damping controls')
    tot_states = 0;
    for i = 1:g.dc.ndci_ud
        ydcimx = g.dc.dci_dc{i,5};
        ydcrmn = g.dc.dci_dc{i,6};
        inv_num = g.dc.dci_dc{i,2};
        st_state = tot_states+1;
        dci_states = g.dc.dci_dc{i,7};
        tot_states = tot_states+dci_states;
        [g.dc.dci_dsig(inv_num,1),g.dc.xdci_dc(st_state:tot_states,1),g.dc.dxdci_dc(st_state:tot_states,1)] =...
            dci_sud(i,1,flag,g.dc.dci_dc{i,1},g.dc.dcid_sig(i,1),ydcimx,ydcimn);
    end
    clear i
end

%% initialize load modulation control
%if ~isempty(lmod_con) % original line - thad
% if statement redundant - used in script... - thad 06/08/20
if ~isempty(g.lmod.lmod_con)
    disp('load modulation')
    lmod(0,1,flag); % removed bus - thad
end
if ~isempty(g.rlmod.rlmod_con)
    disp('reactive load modulation')
    rlmod(0,1,flag); % removed bus - thad
end

%% initialize power modulation control - copied from v2.3 06/01/20 -thad
% Seems like this should be put in a seperate script - thad 06/08/20
if g.pwr.n_pwrmod~=0
    disp('power modulation')
    pwrmod_p(0,1,g.bus.bus,flag);
    pwrmod_q(0,1,g.bus.bus,flag);
    
    [~,~,~,~,Pini,Qini] = pwrmod_dyn([],[],g.bus.bus,g.sys.t,0,0,g.pwr.n_pwrmod);
    
    if (~iscell(Pini) || ~iscell(Qini))
        error('Error in pwrmod_dyn, P_statesIni and P_statesIni must be cells');
    end
    if (any(size(Pini)-[g.pwr.n_pwrmod 1]) || any(size(Qini)-[g.pwr.n_pwrmod 1]))
        error('Dimension error in pwrmod_dyn');
    end
    
    % converted to globals - thad 08/05/20
    g.pwr.pwrmod_p_sigst = cell(g.pwr.n_pwrmod,1);
    g.pwr.pwrmod_q_sigst = g.pwr.pwrmod_p_sigst;
    g.pwr.dpwrmod_p_sigst = g.pwr.pwrmod_p_sigst;
    g.pwr.dpwrmod_q_sigst = g.pwr.pwrmod_p_sigst;
    
    k = size(g.sys.t,2); % length of logged values
    
    for index=1:g.pwr.n_pwrmod
        if ((size(Pini{index},2)~=1) || (size(Qini{index},2)~=1))
            error('Dimension error in pwrmod_dyn');
        end
        % converted to globals - thad 08/05/20
        g.pwr.dpwrmod_p_sigst{index} = zeros(size(Pini{index},1),k);
        g.pwr.pwrmod_p_sigst{index} = Pini{index}*ones(1,k);
        g.pwr.dpwrmod_q_sigst{index} = zeros(size(Qini{index},1),k);
        g.pwr.pwrmod_q_sigst{index} = Qini{index}*ones(1,k);
    end
    clear index dp dq Pini Qini
end

%% initialize non-linear loads
% if statement redundant - used in script... - thad 06/08/20
if ~isempty(g.ncl.load_con)
    disp('non-linear loads')
    vnc = nc_load(g.bus.bus,flag,g.int.Y_ncprf,g.int.Y_ncgprf); % return not used? - thad 07/17/20
else
    g.ncl.nload = 0;
end

%% Initialize AGC
if g.agc.n_agc ~= 0
    calcAreaVals(1,1); % requried for interchange values
    agc(1,0); % initialize
end

%% DC Stuff ? (5/22/20)
if ~isempty(g.dc.dcsp_con)
    % Seems like this should be put in a seperate script - thad 06/08/20
    disp('dc converter specification')
    
    g.bus.bus_sim = g.bus.bus;
    g.bus.bus_int = g.bus.bus_intprf;
    Y1 = g.int.Y_gprf;
    Y2 = g.int.Y_gncprf;
    Y3 = g.int.Y_ncgprf;
    Y4 = g.int.Y_ncprf;
    Vr1 = g.int.V_rgprf;
    Vr2 = g.int.V_rncprf;
    bo = g.int.boprf;
    
    g.k.h_sol = i_simu(1,1,g.k.k_inc, g.k.h, g.bus.bus_sim,Y1,Y2,Y3,Y4,Vr1,Vr2,bo);
    % reinitialize dc controls
    mdc_sig(1);
    dc_cont(0,1,1,g.bus.bus,flag);
    % initialize dc line
    dc_line(0,1,1,g.bus.bus,flag);
end

end%end initNLsim
