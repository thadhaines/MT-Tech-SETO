% s_simu.m
% 9:59 AM 14 June 1999
% An m.file to simulate power system transients
% using the Matlab Power System Toolbox
% This m-file takes the dynamic and load flow data and
% calculates the response of the power system to a fault
% which is specified in a switching file
% see one of the supplied data files (data2a.m) for the
% switching file format

% version 1.9
% Author: Thad Haines
% June 01, 2020
% Modified from PSTv3 to accomodate batch style running
% pwrmod, ivmmod, and gen trip sections from 2.3 copied - untested

% version 1.8
% Author: Rui
% July 19, 2017
% add code for multiple HVDC systems

% version 1.7
% Author: Graham Rogers
% Date September 1999
% Add smple exciter with pi avr - smppi
% version 1.6
% Author: Graham Rogers
% Date: June 1999
% Modification: add user defined hvdc control
%               modify dc so that it is integrated with a sub-multiple time constant

% version 1.5
% Author: Graham Rogers
% Date: December 1998/ January 1999
% Modification: Add svc and tcsc user defined damping controls, and add tcsc model

% version 1.4
% Author: Graham Rogers
% Date:  July 1998
% Modification: Add deltaP/omega filter

% version 1.3
% Author: Graham Rogers
% Date:  June 1998
% Modification: Add hydraulic turbine/governor

% Version: 1.2
% Author:  Graham Rogers
% Date:    August 1997
% Modification : add induction generator
% Version: 1.1
% Author:  Graham Rogers
% Date:    August 1997
% Modification : add modulation for load, exciter and governor
% Version  1.0
% Author:  Graham Rogers
% Date: February 1997
% (c) Copyright: Joe Chow/ Cherry Tree Scientific Software 1991 to 1997 - All rights reserved
%


%%
%clear all
%clear global
% the above clears were removed to allow for running w/o running DataFile.m 5/20/20
% assumes required arrays created before this script runs and DataFile is delted

warning('*** s_simu_Batch Start')

close % close graphics windows
tic % set timer
% plot_now=0;
jay = sqrt(-1);

%% Contents of pst_var copied into this section so that globals highlight - thad
%pst_var % set up global variables (very many)

warning('*** Declare Global Variables')

% Old globals from pst_var

% %     %% system variables - 13
% %     global  basmva basrad syn_ref mach_ref sys_freq
% %     global  bus_v bus_ang psi_re psi_im cur_re cur_im bus_int
% %     % lmon_con not used in non-linear sim...
% %     global  lmon_con

%     %% synchronous machine variables  - 47
%     global  mac_con mac_pot mac_int ibus_con
%     global  mac_ang mac_spd eqprime edprime psikd psikq
%     global  curd curq curdg curqg fldcur
%     global  psidpp psiqpp vex eterm theta ed eq
%     global  pmech pelect qelect
%     global  dmac_ang dmac_spd deqprime dedprime dpsikd dpsikq
%     global  n_mac n_em n_tra n_sub n_ib
%     global  mac_em_idx mac_tra_idx mac_sub_idx mac_ib_idx not_ib_idx
%     global  mac_ib_em mac_ib_tra mac_ib_sub n_ib_em n_ib_tra n_ib_sub
%     global pm_sig n_pm 

    %% excitation system variables - 63
%     global  exc_con exc_pot n_exc
%     global  Efd V_R V_A V_As R_f V_FB V_TR V_B
%     global  dEfd dV_R dV_As dR_f dV_TR
%     global  exc_sig 
%     global smp_idx n_smp dc_idx n_dc  dc2_idx n_dc2 st3_idx n_st3;
%     global smppi_idx n_smppi smppi_TR smppi_TR_idx smppi_no_TR_idx ;
%     global smp_TA smp_TA_idx smp_noTA_idx smp_TB smp_TB_idx smp_noTB_idx;
%     global smp_TR smp_TR_idx smp_no_TR_idx ;
%     global dc_TA dc_TA_idx dc_noTR_idx dc_TB dc_TB_idx dc_noTB_idx;
%     global dc_TE  dc_TE_idx dc_noTE_idx;
%     global dc_TF dc_TF_idx dc_TR dc_TR_idx
%     global st3_TA st3_TA_idx st3_noTA_idx st3_TB st3_TB_idx st3_noTB_idx;
%     global st3_TR st3_TR_idx st3_noTR_idx;
    %% load modulation variables - 7
    %global lmod_con % defined by user
    %global n_lmod lmod_idx % initialized and created in lm_indx
    %global lmod_sig lmod_st dlmod_st % initialized in s_simu
    %global lmod_pot % created/initialized in lmod.m
    % g.lmod.lmod_pot(:,1) = max, g.lmod.lmod_pot(:,2) = min
    %global lmod_data % added by Trudnowski - doesn't appear to be used? maybe in new models?

    % reactive load modulation variables - 7
    %global  rlmod_con n_rlmod rlmod_idx
    %global  rlmod_pot rlmod_st drlmod_st
    %global  rlmod_sig
    
%     %% power injection variables - 10 - g.pwr.
%     global  pwrmod_con n_pwrmod pwrmod_idx
%     global  pwrmod_p_st dpwrmod_p_st
%     global  pwrmod_q_st dpwrmod_q_st
%     global  pwrmod_p_sig pwrmod_q_sig
%     global  pwrmod_data

%     %% non-conforming load variables - 3
%     global  load_con load_pot nload

    %% turbine-governor variables -17
    %global  tg_con tg_pot
    %global  tg1 tg2 tg3 tg4 tg5 dtg1 dtg2 dtg3 dtg4 dtg5
    %global  tg_idx  n_tg tg_sig tgh_idx n_tgh
    
    %% simulation control - 1
%     global sw_con  %scr_con not used

%     %% pss variables - 21
%     global  pss_con pss_pot pss_mb_idx pss_exc_idx
%     global  pss1 pss2 pss3 dpss1 dpss2 dpss3 pss_out
%     global  pss_idx n_pss pss_sp_idx pss_p_idx;
%     global  pss_T  pss_T2 pss_T4 pss_T4_idx  pss_noT4_idx;
% 
%     %% svc variables - 13
%     global  svc_con n_svc svc_idx svc_pot svcll_idx
%     global  svc_sig
%     % svc user defined damping controls
%     global n_dcud dcud_idx svc_dsig
%     global svc_dc % user damping controls?
%     global dxsvc_dc xsvc_dc
%     %states
%     global B_cv B_con
%     %dstates
%     global dB_cv dB_con
% 
%     %% tcsc variables - 10
%     global  tcsc_con n_tcsc tcsvf_idx tcsct_idx
%     global  B_tcsc dB_tcsc
%     global  tcsc_sig tcsc_dsig
%     global  n_tcscud dtcscud_idx  %user defined damping controls
% 	% previous non-globals added as they seem to relavant
% 	global xtcsc_dc dxtcsc_dc td_sig tcscf_idx 
%     global tcsc_dc
%
%     %% induction genertaor variables - 19
%     global  tmig  pig qig vdig vqig  idig iqig igen_con igen_pot
%     global  igen_int igbus n_ig
%     %states
%     global  vdpig vqpig slig
%     %dstates
%     global dvdpig dvqpig dslig
%     % added globals
%     global s_igen
    
%     %% induction motor variables - 21
%     global  tload t_init p_mot q_mot vdmot vqmot  idmot iqmot ind_con ind_pot
%     global  motbus ind_int mld_con n_mot t_mot
%     % states
%     global  vdp vqp slip
%     % dstates
%     global dvdp dvqp dslip
%     % added globals
%     global s_mot
%
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
      
%% Remaining 'loose' globals

    %% ivm variables - 5
    global n_ivm mac_ivm_idx ivmmod_data ivmmod_d_sig ivmmod_e_sig
    
    %% DeltaP/omega filter variables - 21
    global  dpw_con dpw_out dpw_pot dpw_pss_idx dpw_mb_idx dpw_idx n_dpw dpw_Td_idx dpw_Tz_idx
    global  sdpw1 sdpw2 sdpw3 sdpw4 sdpw5 sdpw6
    global  dsdpw1 dsdpw2 dsdpw3 dsdpw4 dsdpw5 dsdpw6

    %% pss design - 3
    global ibus_con  netg_con  stab_con
    
    %% global structured array
    global g


%% meant to be dc globals? 06/08/20 - thad
% look like user defined damping controls that were never implemented
% and intentionally removed/ignored?
g.svc.svc_dc=[];

g.tcsc.tcsc_dc=[];
g.tcsc.n_tcscud = 0;

g.dc.dcr_dc=[];
g.dc.dci_dc=[];

% input data file from m.file
%% 05/20 Edits - thad
% Check for Octave, automatically load compatibility scripe
% Assumes license of Octave will be 'GNU ...'
dataCheck = license;
if all(dataCheck(1:3)=='GNU')
   fprintf('*** Octave detected, loading compatiblity commands and packages...\n')
  octaveComp
else
  fprintf('*** MATLAB detected.\n')
end
clear dataCheck

% account for non-existant DataFile (assumes required arrays created other ways...)
try
   DataFile %Batch name for data file
catch ME
   fprintf('*** Caught ERROR: %s\n',ME.message)
   fprintf('*** Continuing with simulation...\n')
end 
%% run script to handle legacy input to new global g approach
handleNewGlobals

% check for valid dynamic data file
if isempty(g.mac.mac_con)
    error('mac_con is Empty - invalid/incomplete input data.')
end
if isempty(g.sys.sw_con)
    error('sw_con is Empty - simulation has no switching data.')
end

%% Allow Fbase and Sbase to be defined in batch runs
% Handle varaible input system frequency
% assumes fBase defined in DataFile or earlier, sys_freq is defined as global in pst_var.
if ~exist('Fbase','var')
    fprintf('*** Fbase Not defined - assuming 60 Hz base.\n')
    g.sys.sys_freq = 60;
    g.sys.Fbase = 60;
elseif isnumeric(Fbase)
    fprintf('*** Fbase found - Frequency base is set to %3.3f Hz\n', Fbase) 
    g.sys.sys_freq = Fbase;
    g.sys.Fbase = Fbase;
end

% Handle variable input base MVA
% assumes Sbase defined in DataFile or earlier, basmva is defined as global in pst_var.
if ~exist('Sbase','var')
    fprintf('*** Sbase Not defined - assuming 100 MVA base.\n')
    g.sys.basmva = 100;
elseif isnumeric(Sbase)
    fprintf('*** Sbase found - Power base is set to %3.3f MVA\n', Sbase) 
    g.sys.basmva = Sbase;
end

%% other init operations
g.sys.basrad = 2*pi*g.sys.sys_freq; % default system frequency is 60 Hz
g.sys.syn_ref = 0 ;     % synchronous reference frame
g.mac.ibus_con = []; % ignore infinite buses in transient simulation % should be global? -thad 06/15/20

%% solve for loadflow - loadflow parameter
warning('*** Solve initial loadflow')
if isempty(g.dc.dcsp_con)
    % AC power flow
    g.dc.n_conv = 0;
    g.dc.n_dcl = 0;
    g.dc.ndcr_ud=0;
    g.dc.ndci_ud=0;
    tol = 1e-9;   % tolerance for convergence
    iter_max = 30; % maximum number of iterations
    acc = 1.0;   % acceleration factor
    [bus_sol,line,line_flw] = loadflow(bus,line,tol,iter_max,acc,'n',2);
    bus = bus_sol;  % solved loadflow solution needed for initialization
    save sim_fle.mat bus line
else
    % Has HVDC, use DC load flow
    [bus_sol,line,line_flw,rec_par,inv_par, line_par] = lfdcs(bus,line,g.dc.dci_dc,g.dc.dcr_dc);
    bus = bus_sol;
    save sim_fle.mat bus line rec_par  inv_par line_par
end

%% set indexes
% note: dc index set in dc load flow
mac_indx();
exc_indx();
tg_indx(); % functionalized 06/05/20 - thad
dpwf_indx;
pss_indx;
svc_indx(); % assigned svc_dc to global (damping control?)
tcsc_indx();
lm_indx;
rlm_indx();
pwrmod_indx(bus); 

g.ind.n_mot = size(g.ind.ind_con,1); % inductive motors
g.igen.n_ig = size(g.igen.igen_con,1); % inductive generators

if isempty(g.ind.n_mot)
    g.ind.n_mot = 0;
end
if isempty(g.igen.n_ig)
    g.igen.n_ig = 0; 
end

ntot = g.mac.n_mac+g.ind.n_mot+g.igen.n_ig;
ngm = g.mac.n_mac + g.ind.n_mot;

g.mac.n_pm = g.mac.n_mac; % into an init?

%% Make sure bus max/min Q is the same as the pwrmod_con max/min Q - moved to place after pwrmod is counted (never actually executed prior...) -thad 06/30/20
if ~isempty(g.pwr.n_pwrmod)
    for kk=1:g.pwr.n_pwrmod
        n = find(g.pwr.pwrmod_con(kk,1)==bus(:,1));
        bus(n,11:12) = g.pwr.pwrmod_con(kk,6:7);
    end
    clear kk n
end


%% construct simulation switching sequence as defined in sw_con
warning('*** Initialize time and switching variables')
tswitch(1) = g.sys.sw_con(1,1);
k = 1;
kdc = 1;
n_switch = length(g.sys.sw_con(:,1));
k_inc = zeros(n_switch-1,1);
k_incdc = k_inc;
t_switch = zeros(n_switch,1);
h = t_switch;
h_dc = h;

for sw_count = 1:n_switch-1
    h(sw_count) = g.sys.sw_con(sw_count,7);%specified time step
    
    if h(sw_count)==0
        h(sw_count) = 0.01;
    end % default time step
    
    k_inc(sw_count) = fix((g.sys.sw_con(sw_count+1,1)-g.sys.sw_con(sw_count,1))/h(sw_count));%nearest lower integer
    
    if k_inc(sw_count)==0
        k_inc(sw_count)=1;
    end% minimum 1
    
    h(sw_count) = (g.sys.sw_con(sw_count+1,1)-g.sys.sw_con(sw_count,1))/k_inc(sw_count);%step length
    h_dc(sw_count) = h(sw_count)/10;
    k_incdc(sw_count) = 10*k_inc(sw_count);
    t_switch(sw_count+1) = t_switch(sw_count) +  k_inc(sw_count)*h(sw_count);
    t(k:k-1+k_inc(sw_count)) = t_switch(sw_count):h(sw_count):t_switch(sw_count+1)-h(sw_count);
    g.dc.t_dc(kdc:kdc-1+k_incdc(sw_count)) = t_switch(sw_count):h_dc(sw_count):t_switch(sw_count+1)-h_dc(sw_count);
    k = k+k_inc(sw_count);
    kdc = kdc+k_incdc(sw_count);
end

% time for dc - multi-rate...
g.dc.t_dc(kdc)=g.dc.t_dc(kdc-1)+h_dc(sw_count);
for kk=1:10
    kdc=kdc+1;
    g.dc.t_dc(kdc)=g.dc.t_dc(kdc-1)+h_dc(sw_count);
end

k = sum(k_inc)+1; % k is the total number of time steps in the simulation

t(k) = g.sys.sw_con(n_switch,1);
n = size(g.mac.mac_con, 1) ;
n_bus = length(bus(:,1));

% add time array t to global g - thad
g.sys.t = t;

%% create zero matrices for variables to make algorithm more efficient?
warning('*** Initialize zero matricies...')
z = zeros(n,k);
z1 = zeros(1,k);
zm = zeros(1,k);
if g.ind.n_mot>1
    zm = zeros(g.ind.n_mot,k);
end

zig = zeros(1,k);
if g.igen.n_ig>1
    zig = zeros(g.igen.n_ig,k);
end

zdc = zeros(2,kdc);
if g.dc.n_conv>2
    zdc = zeros(g.dc.n_conv,kdc);
end

zdcl = zeros(1,kdc);
if g.dc.n_dcl>1
    zdcl=zeros(g.dc.n_dcl,kdc);
end

%% set dc parameters   (initialize zeros? 5/14/20)
g.dc.Vdc = zeros(g.dc.n_conv,kdc);
g.dc.i_dc = zdc;
g.dc.P_dc = z; 
g.dc.cur_ord = z;
g.dc.alpha = zdcl;
g.dc.gamma = zdcl;
g.dc.dc_sig = zeros(g.dc.n_conv,k);
g.dc.dc_dsig = zeros(g.dc.n_conv,k);
g.dc.dcr_dsig = zeros(g.dc.n_dcl,k);
g.dc.dci_dsig = zeros(g.dc.n_dcl,k);
g.dc.i_dcr = zdcl; 
g.dc.i_dci = zdcl; 
g.dc.v_dcc = zdcl;
g.dc.di_dcr = zdcl; 
g.dc.di_dci = zdcl; 
g.dc.dv_dcc = zdcl;
g.dc.v_conr = zdcl; 
g.dc.v_coni = zdcl; 
g.dc.dv_conr = zdcl; 
g.dc.dv_coni = zdcl;

if g.dc.n_conv~=0
    % Modified by Rui on Oct. 5, 2016
    for ihvdc_count=1:g.dc.n_dcl
        ire = g.dc.r_idx(ihvdc_count);
        g.dc.Vdc(ire,:) = rec_par(ihvdc_count,2);
        g.dc.i_dc(ire,:) = line_par(ihvdc_count);    % for PDCI
        g.dc.i_dcr(ihvdc_count,:) = g.dc.i_dc(ire,:);
        g.dc.alpha(ihvdc_count,:) = rec_par(ihvdc_count,1)*pi/180;
    end
    
    for ihvdc_count=1:g.dc.n_dcl
        iin=g.dc.i_idx(ihvdc_count);
        g.dc.Vdc(iin,:)=inv_par(ihvdc_count,2);
        g.dc.i_dc(iin,:) = line_par(ihvdc_count);
        g.dc.i_dci(ihvdc_count,:) = g.dc.i_dc(iin,:);
        g.dc.gamma(ihvdc_count,:) = inv_par(ihvdc_count,1)*pi/180;
    end
    % end modification by Rui
    
    g.dc.Vdc(g.dc.r_idx,:) = rec_par(:,2); 
    g.dc.Vdc(g.dc.i_idx,:) = inv_par(:,2);
    g.dc.i_dc(g.dc.r_idx,:) = line_par; 
    g.dc.i_dc(g.dc.i_idx,:) = line_par;
    g.dc.i_dcr(:,:) = g.dc.i_dc(g.dc.r_idx,:); 
    g.dc.i_dci(:,:) = g.dc.i_dc(g.dc.i_idx,:);
    g.dc.alpha(:,:) = rec_par(:,1)*pi/180;
    g.dc.gamma(:,:) = inv_par(:,1)*pi/180;
    
    if g.dc.ndcr_ud~=0
        for j = 1:g.dc.ndcr_ud
            sv = get(g.dc.dcr_dc{j,1});
            if j==1
                g.dc.xdcr_dc =zeros(sv.NumStates,kdc); 
                g.dc.dxdcr_dc = zeros(sv.NumStates,kdc);
            else
                g.dc.xdcr_dc = [g.dc.xdcr_dc;zeros(sv.NumStates,kdc)];
                g.dc.dxdcr_dc = [g.dc.dxdcr_dc;zeros(sv.NumStates,kdc)];
            end
        end
        g.dc.dcrd_sig=zeros(g.dc.ndcr_ud,k);
        g.dc.angdcr = zeros(g.dc.ndcr_ud,k);
    else
        g.dc.xdcr_dc = zeros(1,kdc);
        g.dc.dxdcr_dc = zeros(1,kdc);
        g.dc.dcrd_sig = zeros(1,k);
    end
    
    if g.dc.ndci_ud~=0
        for j = 1:g.dc.ndci_ud
            sv = get(g.dc.dci_dc{j,1});
            if j==1
                g.dc.xdci_dc =zeros(sv.NumStates,kdc);
                g.dc.dxdci_dc = zeros(sv.NumStates,kdc);
            else
                g.dc.xdci_dc = [g.svc.xsvc_dc;zeros(sv.NumStates,kdc)];
                g.dc.dxdci_dc = [g.svc.dxsvc_dc;zeros(sv.NumStates,kdc)];
            end
        end
        g.dc.dcid_sig=zeros(g.dc.ndcr_ud,k);
        g.dc.angdci = zeros(ndci_dc,k);
    else
        g.dc.xdci_dc = zeros(1,kdc);
        g.dc.dxdci_dc = zeros(1,kdc);
        g.dc.dcid_sig = zeros(1,k);
    end
else
    g.dc.xdcr_dc = zeros(1,kdc);
    g.dc.dxdcr_dc = zeros(1,kdc);
    g.dc.xdci_dc = zeros(1,kdc);
    g.dc.dxdci_dc = zeros(1,kdc);
end

%% End of DC specific stuff? - continuation of zero intialization - 06/01/20 -thad

v_p = z1;
mac_ref = z1;  % unsure of this use
sys_ref = z1;   % unsure of this use - thad 07/02/20

g.sys.bus_v = zeros(n_bus+1,k);

g.sys.cur_re = z; 
g.sys.cur_im = z; 
g.sys.psi_re = z; 
g.sys.psi_im = z;

g.sys.theta = zeros(n_bus+1,k);

g.mac.mac_ang = z; 
g.mac.mac_spd = z; 
g.mac.dmac_ang = z; 
g.mac.dmac_spd = z;
g.mac.pmech = z; 
g.mac.pelect = z; 
g.mac.edprime = z; 
g.mac.eqprime = z; 
g.mac.dedprime = z; 
g.mac.deqprime = z;
g.mac.psikd = z; 
g.mac.psikq = z; 
g.mac.dpsikd = z; 
g.mac.dpsikq = z;
g.mac.pm_sig = z;
g.mac.curd = z; 
g.mac.curq = z; 
g.mac.curdg = z; 
g.mac.curqg = z; 
g.mac.fldcur = z;
g.mac.ed = z; 
g.mac.eq = z; 
g.mac.eterm = z; 
g.mac.qelect = z;
g.mac.vex = z; 

% modified to use global g - thad 06/05/20
totGovs = g.tg.n_tg + g.tg.n_tgh;
if totGovs ~= 0
    z_tg = zeros(totGovs,k);
else
    z_tg = zeros(1,k);
end
% seems unrequired to create these zeros if totGovs == 0... -thad
g.tg.tg1 = z_tg; 
g.tg.tg2 = z_tg; 
g.tg.tg3 = z_tg; 
g.tg.tg4 = z_tg; 
g.tg.tg5 = z_tg;
g.tg.dtg1 = z_tg; 
g.tg.dtg2 = z_tg; 
g.tg.dtg3 = z_tg;
g.tg.dtg4 = z_tg; 
g.tg.dtg5 = z_tg;
g.tg.tg_sig = z_tg;
clear totGovs z_tg
%

z_pss = zeros(1,k);
if g.pss.n_pss~=0
    z_pss = zeros(g.pss.n_pss,k);
end

g.pss.pss1 = z_pss; 
g.pss.pss2 = z_pss; 
g.pss.pss3 = z_pss;
g.pss.dpss1 = z_pss; 
g.pss.dpss2 = z_pss; 
g.pss.dpss3 = z_pss;

z_dpw = zeros(1,k);
if n_dpw~=0
    z_dpw = zeros(n_dpw,k);
end

sdpw1 = z_dpw; 
sdpw2 = z_dpw; 
sdpw3 = z_dpw; 
sdpw4 = z_dpw; 
sdpw5 = z_dpw; 
sdpw6 = z_dpw; 
dpw_out = z_dpw;
dsdpw1 = z_dpw; 
dsdpw2 = z_dpw; 
dsdpw3 = z_dpw; 
dsdpw4 = z_dpw; 
dsdpw5 = z_dpw; 
dsdpw6 = z_dpw;


%% exciter zero init
ze = zeros(1,k);
if g.exc.n_exc~=0
    ze = zeros(g.exc.n_exc,k);
end

g.exc.V_B = ze;
g.exc.exc_sig = ze;
g.exc.V_TR = ze; 
g.exc.V_R = ze; 
g.exc.V_A = ze; 
g.exc.V_As = ze; 
g.exc.Efd = ze; 
g.exc.R_f = ze;
g.exc.dV_TR = ze; 
g.exc.dV_R = ze; 
g.exc.dV_As = ze; 
g.exc.dEfd = ze; 
g.exc.dR_f = ze;

g.pss.pss_out = ze; % seems related to pss more than exciters - thad 06/17/20

%% inductive motor zeros
g.ind.vdp = zm; 
g.ind.vqp = zm; 
g.ind.slip = zm;
g.ind.dvdp = zm; 
g.ind.dvqp = zm; 
g.ind.dslip = zm;
g.ind.s_mot = zm; % added to global 07/13/20 - thad
g.ind.p_mot = zm; 
g.ind.q_mot = zm;

g.igen.vdpig = zig; 
g.igen.vqpig = zig; 
g.igen.slig = zig;
g.igen.dvdpig = zig; 
g.igen.dvqpig = zig; 
g.igen.dslig = zig;
g.igen.s_igen = zig; % changed to global -thad 07/13/20
g.igen.pig = zig; 
g.igen.qig = zig; 
g.igen.tmig = zig;

if g.svc.n_svc~=0
    svcZeros = zeros(g.svc.n_svc,k); 
    g.svc.B_cv = svcZeros;
    g.svc.dB_cv = svcZeros;
    g.svc.svc_sig = svcZeros;
    g.svc.svc_dsig= svcZeros;
    g.svc.B_con = svcZeros;
    g.svc.dB_con= svcZeros;
    
    % Unsure of this - DC svc? -thad 07/08/20
    % damping control user defined?
    if g.svc.n_dcud~=0
        d_sig = zeros(g.svc.n_dcud,k);
        for j = 1:g.svc.n_dcud
            sv = get(g.svc.svc_dc{j,1});
            if j==1
                g.svc.xsvc_dc =zeros(sv.NumStates,k);
                g.svc.dxsvc_dc = zeros(sv.NumStates,k);
            else
                g.svc.xsvc_dc = [g.svc.xsvc_dc;zeros(sv.NumStates,k)];
                g.svc.dxsvc_dc = [g.svc.dxsvc_dc;zeros(sv.NumStates,k)];
            end
        end
        
    else
        g.svc.xsvc_dc = zeros(1,k);
        g.svc.dxsvc_dc = zeros(1,k);
    end
else
    g.svc.B_cv = zeros(1,k);
    g.svc.dB_cv = zeros(1,k); 
    g.svc.svc_sig = zeros(1,k);
    g.svc.svc_dsig = zeros(1,k);
    g.svc.B_con = zeros(1,k);
    g.svc.dB_con=zeros(1,k);
    
    % DC SVC?
    g.svc.xsvc_dc = zeros(1,k);   % supposed to be global? - thad 06/03/20
    g.svc.dxsvc_dc = zeros(1,k);  % supposed to be global? - thad 06/03/20
    d_sig = zeros(1,k);     % supposed to be global? - thad 06/03/20
end

if g.tcsc.n_tcsc~=0
    g.tcsc.B_tcsc = zeros(g.tcsc.n_tcsc,k); 
    g.tcsc.dB_tcsc = zeros(g.tcsc.n_tcsc,k);
    g.tcsc.tcsc_sig = zeros(g.tcsc.n_tcsc,k);
    g.tcsc.tcsc_dsig=zeros(g.tcsc.n_tcsc,k);
    
    if g.tcsc.n_tcscud~=0
        g.tcsc.td_sig = zeros(g.tcsc.n_tcscud,k);%input to tcsc damping control
        for j = 1:g.tcsc.n_tcscud
            sv = get(g.tcsc.tcsc_dc{j,1});% damping control state space object
            if j==1
                g.tcsc.xtcsc_dc =zeros(sv.NumStates,k); % tcsc damping control states
                g.tcsc.dxtcsc_dc = zeros(sv.NumStates,k);% tcsc dc rates of chage of states
            else
                g.tcsc.xtcsc_dc = [g.tcsc.xtcsc_dc;zeros(sv.NumStates,k)];% in order of damping controls
                g.tcsc.dxtcsc_dc = [g.tcsc.dxtcsc_dc;zeros(sv.NumStates,k)];
            end
        end
    else
        g.tcsc.xtcsc_dc = zeros(1,k);
        g.tcsc.dxtcsc_dc = zeros(1,k);
    end
else
    g.tcsc.B_tcsc = zeros(1,k);
    g.tcsc.dB_tcsc = zeros(1,k); 
    g.tcsc.tcsc_sig = zeros(1,k);
    g.tcsc.tcsc_dsig = zeros(1,k);
    g.tcsc.xtcsc_dc = zeros(1,k);
    g.tcsc.dxtcsc_dc = zeros(1,k);
    g.tcsc.td_sig = zeros(1,k);
end


if g.lmod.n_lmod ~= 0 
    % initialize zeros for all states and signals associated with lmod
    g.lmod.lmod_st = zeros(g.lmod.n_lmod,k);
    g.lmod.dlmod_st = g.lmod.lmod_st;
    g.lmod.lmod_sig = g.lmod.lmod_st;
else
    % initialize single row of zeros ( may be un necessary) - thad
    g.lmod.lmod_st = zeros(1,k);
    g.lmod.dlmod_st = g.lmod.lmod_st;
    g.lmod.lmod_sig = g.lmod.lmod_st;
end

if g.rlmod.n_rlmod ~= 0
    g.rlmod.rlmod_st = zeros(g.rlmod.n_rlmod,k); 
    g.rlmod.drlmod_st = g.rlmod.rlmod_st; 
    g.rlmod.rlmod_sig = g.rlmod.rlmod_st;
else
    g.rlmod.rlmod_st = zeros(1,k); 
    g.rlmod.drlmod_st = g.rlmod.rlmod_st; 
    g.rlmod.rlmod_sig = g.rlmod.rlmod_st;
end

%% Initialize pwrmod
if g.pwr.n_pwrmod ~= 0
    g.pwr.pwrmod_p_st = zeros(g.pwr.n_pwrmod,k);
    g.pwr.dpwrmod_p_st = g.pwr.pwrmod_p_st;
    g.pwr.pwrmod_p_sig = g.pwr.pwrmod_p_st;
    g.pwr.pwrmod_q_st = zeros(g.pwr.n_pwrmod,k);
    g.pwr.dpwrmod_q_st = g.pwr.pwrmod_q_st;
    g.pwr.pwrmod_q_sig = g.pwr.pwrmod_q_st;
else
    g.pwr.pwrmod_p_st = zeros(1,k);
    g.pwr.dpwrmod_p_st = g.pwr.pwrmod_p_st;
    g.pwr.pwrmod_p_sig = g.pwr.pwrmod_p_st;
    g.pwr.pwrmod_q_st = zeros(1,k);
    g.pwr.dpwrmod_q_st = g.pwr.pwrmod_q_st;
    g.pwr.pwrmod_q_sig = g.pwr.pwrmod_q_st;
end

%% Initialize ivmmod sigs
if n_ivm ~= 0
    ivmmod_d_sig = zeros(n_ivm,k);
    ivmmod_e_sig = zeros(n_ivm,k);
else
    ivmmod_d_sig = zeros(1,k);
    ivmmod_e_sig = zeros(1,k);
end

g.sys.sys_freq = ones(1,k); % replaces variable for base frequency input... 5/21/20

%% step 1: construct reduced Y matrices
warning('*** Initialize Y matrix (matracies?) and Dynamic Models')
disp('constructing reduced y matrices')
disp('initializing motor,induction generator, svc and dc control models')

bus = mac_ind(0,1,bus,0);% initialize induction motor
bus = mac_igen(0,1,bus,0); % initialize induction generator
bus = svc(0,1,bus,0);%initialize svc
dc_cont(0,1,1,bus,0);% initialize dc controls
% this has to be done before red_ybus is used since the motor and svc
% initialization alters the bus matrix and dc parameters are required

y_switch % calculates the reduced y matrices for the different switching conditions

disp('initializing other models...')

%% step 2: initialization
g.sys.theta(1:n_bus,1) = bus(:,3)*pi/180;
g.sys.bus_v(1:n_bus,1) = bus(:,2).*exp(jay*g.sys.theta(1:n_bus,1));

if g.svc.n_dcud ~=0 % Seems like this should be put in a seperate script - thad 06/08/20
    %% calculate the initial magnitude of line current for svc damping controls
    for j=1:g.svc.n_dcud
        l_num = g.svc.svc_dc{j,3};
        svc_num = g.svc.svc_dc{j,2};
        from_bus = g.sys.bus_int(line(l_num,1)); 
        to_bus = g.sys.bus_int(line(l_num,2));
        svc_bn = g.sys.bus_int(g.svc.svc_con(svc_num,2));
        
        if svc_bn~= from_bus&& svc_bn  ~= to_bus
            error('the svc is not at the end of the specified line');
        end
        
        V1 = g.sys.bus_v(from_bus,1);
        V2 = g.sys.bus_v(to_bus,1);
        R = line(l_num,3);
        X = line(l_num,4);
        B = line(l_num,5);
        g.dc.tap = line(l_num,6);
        phi = line(l_num,7);
        [l_if,l_it] = line_cur(V1,V2,R,X,B,g.dc.tap,phi);
        l_if0(j)=l_if;
        l_it0(j)=l_it;
        
        if svc_bn == from_bus
            d_sig(j,1)=abs(l_if);
        elseif svc_bn==to_bus
            d_sig(j,1)=abs(l_it);
        end
    end
end

if g.tcsc.n_tcscud ~=0 % Seems like this should be put in a seperate script - thad 06/08/20
    %% calculate the initial magnitude of bus voltage magnitude for tcsc damping controls
    for j=1:g.tcsc.n_tcscud
        b_num = g.tcsc.tcsc_dc{j,3};
        tcsc_num = g.tcsc.tcsc_dc{j,2};
        g.tcsc.td_sig(j,1) =abs (g.sys.bus_v(g.sys.bus_int(b_num),1));
    end
end

if g.dc.n_conv~=0 % Seems like this should be put in a seperate script - thad 06/08/20
    %% change dc buses from LT to HT
    Pr = bus(g.dc.rec_ac_bus,6);
    Pi = bus(g.dc.inv_ac_bus,6);
    Qr = bus(g.dc.rec_ac_bus,7);
    Qi = bus(g.dc.inv_ac_bus,7);
    VLT= g.sys.bus_v(g.dc.ac_bus,1);
    i_acr = (Pr-jay*Qr)./conj(VLT(g.dc.r_idx));
    i_aci = (Pi - jay*Qi)./conj(VLT(g.dc.i_idx));
    IHT(g.dc.r_idx,1)=i_acr;
    IHT(g.dc.i_idx,1)=i_aci;
    g.dc.VHT(g.dc.r_idx,1) = (VLT(g.dc.r_idx) + jay*g.dc.dcc_pot(:,2).*i_acr);
    g.dc.VHT(g.dc.i_idx,1) = (VLT(g.dc.i_idx) + jay*g.dc.dcc_pot(:,4).*i_aci);
    g.sys.bus_v(g.dc.ac_bus,1) = g.dc.VHT;
    g.sys.theta(g.dc.ac_bus,1) = angle(g.sys.bus_v(g.dc.ac_bus,1));
    % modify the bus matrix to the HT buses
    bus(g.dc.ac_bus,2) = abs(g.sys.bus_v(g.dc.ac_bus,1));
    bus(g.dc.ac_bus,3) = g.sys.theta(g.dc.ac_bus,1)*180/pi;
    SHT = g.dc.VHT.*conj(IHT);
    bus(g.dc.ac_bus,6) = real(SHT);
    bus(g.dc.ac_bus,7) = imag(SHT);
    
    if g.dc.ndcr_ud~=0 % Seems like this should be put in a seperate script - thad 06/08/20
        % calculate the initial value of bus angles rectifier user defined control
        for j = 1:g.dc.ndcr_ud
            b_num1 = g.dc.dcr_dc{j,3};
            b_num2 = g.dc.dcr_dc{j,4};
            conv_num = g.dc.dcr_dc{j,2};
            g.dc.angdcr(j,:) = g.sys.theta(g.sys.bus_int(b_num1),1)-g.sys.theta(g.sys.bus_int(b_num2),1);
            g.dc.dcrd_sig(j,:)=g.dc.angdcr(j,:);
        end
    end
    if g.dc.ndci_ud~=0 % Seems like this should be put in a seperate script - thad 06/08/20
        % calculate the initial value of bus angles inverter user defined control
        for j = 1:g.dc.ndci_ud
            b_num1 = g.dc.dci_dc{j,3};
            b_num2 = g.dc.dci_dc{j,4};
            conv_num = g.dc.dci_dc{j,2};
            g.dc.angdci(j,:) = g.sys.theta(g.sys.bus_int(b_num1),1)-g.sys.theta(g.sys.bus_int(b_num2),1);
            g.dc.dcid_sig(j,:) = g.dc.angdci(j,:);
        end
    end
end

%% Flag = 0 == Initialization
warning('*** Dynamic model initialization via functions/scripts:')
flag = 0;
g.sys.bus_int = bus_intprf;% pre-fault system

disp('generators')
mac_sub(0,1,bus,flag); % first 
mac_tra(0,1,bus,flag);
mac_em(0,1,bus,flag);
mac_ivm(0,1,bus,flag); % ivm - thad 06/01/20

disp('generator controls')
dpwf(0,1,flag);
pss(0,1,flag);

% exciters
smpexc(0,1,flag);
smppi(0,1,flag);
exc_st3(0,1,flag);
exc_dc12(0,1,flag);

tg(0,1,flag); % modified 06/05/20 to global g

tg_hydro(0,1,bus,flag);

%% initialize ivm modulation control - added from v2.3 06/01/20 - thad
% Seems like this should be put in a seperate script - thad 06/08/20
if n_ivm~=0
    disp('ivm modulation')
    [~,~,~,~,Dini,Eini] = ivmmod_dyn([],[],bus,t,1,flag);
    
    if (~iscell(Dini) || ~iscell(Eini))
        error('Error in ivmmod_dyn, initial states must be cells'); 
    end
    if (any(size(Dini)-[n_ivm 1]) || any(size(Eini)-[n_ivm 1]))
        error('Dimension error in ivmmod_dyn'); 
    end
    
    ivmmod_d_sigst = cell(n_ivm,1);
    ivmmod_e_sigst = ivmmod_d_sigst;
    divmmod_d_sigst = ivmmod_d_sigst;
    divmmod_e_sigst = ivmmod_d_sigst;
    
    for index=1:n_ivm
        if ((size(Dini{index},2)~=1) || (size(Eini{index},2)~=1))
            error('Dimension error in ivmmod_dyn'); 
        end
        divmmod_d_sigst{index} = zeros(size(Dini{index},1),k);
        ivmmod_d_sigst{index} = Dini{index}*ones(1,k);
        divmmod_e_sigst{index} = zeros(size(Eini{index},1),k);
        ivmmod_e_sigst{index} = Eini{index}*ones(1,k);
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
    pwrmod_p(0,1,bus,flag);
    pwrmod_q(0,1,bus,flag);
    [~,~,~,~,Pini,Qini] = pwrmod_dyn([],[],bus,t,0,0,g.pwr.n_pwrmod);
    if (~iscell(Pini) || ~iscell(Qini))
        error('Error in pwrmod_dyn, P_statesIni and P_statesIni must be cells');
    end
    if (any(size(Pini)-[g.pwr.n_pwrmod 1]) || any(size(Qini)-[g.pwr.n_pwrmod 1]))
        error('Dimension error in pwrmod_dyn'); 
    end
    pwrmod_p_sigst = cell(g.pwr.n_pwrmod,1);
    pwrmod_q_sigst = pwrmod_p_sigst;
    dpwrmod_p_sigst = pwrmod_p_sigst;
    dpwrmod_q_sigst = pwrmod_p_sigst;
    for index=1:g.pwr.n_pwrmod
        if ((size(Pini{index},2)~=1) || (size(Qini{index},2)~=1))
            error('Dimension error in pwrmod_dyn'); 
        end
        dpwrmod_p_sigst{index} = zeros(size(Pini{index},1),k); 
        pwrmod_p_sigst{index} = Pini{index}*ones(1,k);
        dpwrmod_q_sigst{index} = zeros(size(Qini{index},1),k); 
        pwrmod_q_sigst{index} = Qini{index}*ones(1,k);
    end
    clear index dp dq Pini Qini
end

%% initialize non-linear loads
% if statement redundant - used in script... - thad 06/08/20
if ~isempty(g.ncl.load_con)
    disp('non-linear loads')
    vnc = nc_load(bus,flag,Y_ncprf,Y_ncgprf);
else
    g.ncl.nload = 0;
end

%% DC Stuff ? (5/22/20)
if ~isempty(g.dc.dcsp_con)
% Seems like this should be put in a seperate script - thad 06/08/20
    disp('dc converter specification')
    
    bus_sim = bus;
    g.sys.bus_int = bus_intprf;
    Y1 = Y_gprf;
    Y2 = Y_gncprf;
    Y3 = Y_ncgprf;
    Y4 = Y_ncprf;
    Vr1 = V_rgprf;
    Vr2 = V_rncprf;
    bo = boprf;
    
    h_sol = i_simu(1,1,k_inc,h,bus_sim,Y1,Y2,Y3,Y4,Vr1,Vr2,bo);
    % reinitialize dc controls
    mdc_sig(1);
    dc_cont(0,1,1,bus,flag);
    % initialize dc line
    dc_line(0,1,1,bus,flag);
end

H_sum = sum(g.mac.mac_con(:,16)./g.mac.mac_pot(:,1));
%% step 3: perform a predictor-corrector integration
%
% This doesn't seem to be very predictor correctory.... 5/15/20
% more to do with time/data point increment

kt = 0;
ks = 1;

k_tot = sum(k_inc);
lswitch = length(k_inc);
ktmax = k_tot-k_inc(lswitch);
bus_sim = bus;

% added from v2.3 06/01/20 - thad
mac_trip_flags = false(g.mac.n_mac,1);
mac_trip_states = 0;

%% Simulation loop start
warning('*** Simulation Loop Start')
while (kt<=ktmax)
    k_start = kt+1;
    
    if kt==ktmax
        k_end = kt + k_inc(ks);
    else
        k_end = kt + k_inc(ks) + 1;
    end
    
    for k = k_start:k_end
        %% step 3a: network solution
        
        % display k and t at k_inc and every ...th step - thad
        if ( mod(k,50)==0 ) || k == 1 || k == k_end
            fprintf('*** k = %5d, \tt(k) = %7.4f\n',k,t(k)) % DEBUG
        end
        
        % mach_ref(k) = mac_ang(syn_ref,k);
        g.sys.mach_ref(k) = 0;
        g.mac.pmech(:,k+1) = g.mac.pmech(:,k);
        g.igen.tmig(:,k+1) = g.igen.tmig(:,k);
        
        if g.dc.n_conv~=0
            g.dc.cur_ord(:,k+1) = g.dc.cur_ord(:,k);
        end
        
        % Trip gen - Copied from v2.3 06/01/20 - thad
        [f,mac_trip_states] = mac_trip_logic(mac_trip_flags,mac_trip_states,t,k);
        mac_trip_flags = mac_trip_flags | f;
        
        %% Flag = 1    
        flag = 1;
        timestep = int2str(k); % not used? 06/09/20
        % network-machine interface
        mac_ind(0,k,bus_sim,flag);
        mac_igen(0,k,bus_sim,flag);
        mac_sub(0,k,bus_sim,flag);
        mac_tra(0,k,bus_sim,flag);
        mac_em(0,k,bus_sim,flag);
        mac_ivm(0,k,bus_sim,flag); 
        
        mdc_sig(k); % dc controls mod signals
        dc_cont(0,k,10*(k-1)+1,bus_sim,flag); % Models the action of HVDC link pole controllers
        
        %% Calculate current injections and bus voltages and angles
        if k >= sum(k_inc(1:3))+1
            %% fault cleared
            line_sim = line_pf2;
            bus_sim = bus_pf2;
            g.sys.bus_int = bus_intpf2;
            Y1 = Y_gpf2;
            Y2 = Y_gncpf2;
            Y3 = Y_ncgpf2;
            Y4 = Y_ncpf2;
            Vr1 = V_rgpf2;
            Vr2 = V_rncpf2;
            bo = bopf2;
            % i_simu forms the network interface variables
            %h_sol = i_simu(k,ks,k_inc,h,bus_sim,Y1,Y2,Y3,Y4,Vr1,Vr2,bo);
            % duplicate call?
            % h_sol calculated after this 'if' block...
            
        elseif k >=sum(k_inc(1:2))+1
            %% near bus cleared
            line_sim = line_pf1;
            bus_sim = bus_pf1;
            g.sys.bus_int = bus_intpf1;
            Y1 = Y_gpf1;
            Y2 = Y_gncpf1;
            Y3 = Y_ncgpf1;
            Y4 = Y_ncpf1;
            Vr1 = V_rgpf1;
            Vr2 = V_rncpf1;
            bo = bopf1;
            
            %h_sol = i_simu(k,ks,k_inc,h,bus_sim,Y1,Y2,Y3,Y4,Vr1,Vr2,bo);
            
        elseif k>=k_inc(1)+1
            %% fault applied
            line_sim = line_f;
            bus_sim = bus_f;
            g.sys.bus_int = bus_intf;
            Y1 = Y_gf;
            Y2 = Y_gncf;
            Y3 = Y_ncgf;
            Y4 = Y_ncf;
            Vr1 = V_rgf;
            Vr2 = V_rncf;
            bo = bof;
            
            %h_sol = i_simu(k,ks,k_inc,h,bus_sim,Y1,Y2,Y3,Y4,Vr1,Vr2,bo);
            
        elseif k<k_inc(1)+1
            %% pre fault
            line_sim = line;
            bus_sim = bus;
            g.sys.bus_int = bus_intprf;
            Y1 = Y_gprf;
            Y2 = Y_gncprf;
            Y3 = Y_ncgprf;
            Y4 = Y_ncprf;
            Vr1 = V_rgprf;
            Vr2 = V_rncprf;
            bo = boprf;
            
            %h_sol = i_simu(k,ks,k_inc,h,bus_sim,Y1,Y2,Y3,Y4,Vr1,Vr2,bo);
        end
        
        %% apply gen trip - added from v2.3 - 06/01/20 - thad
        if sum(mac_trip_flags)>0.5
            genBuses = g.mac.mac_con(mac_trip_flags==1,2);
            for kB=1:length(genBuses)
                nL = find(genBuses(kB)==line_sim(:,1) | genBuses(kB)==line_sim(:,2));
                if isempty(nL); error(' '); end
                line_sim(nL,4) = 1e7; %make reactance infinity
            end
            [Y1,Y2,Y3,Y4,Vr1,Vr2,bo] = red_ybus(bus_sim,line_sim);
            clear nL kB genBuses
        end
        
        %% solve
        h_sol = i_simu(k,ks,k_inc,h,bus_sim,Y1,Y2,Y3,Y4,Vr1,Vr2,bo);
        
        %% HVDC
        if g.dc.ndcr_ud~=0
            % calculate the new value of bus angles rectifier user defined control
            tot_states=0;
            for jj = 1:g.dc.ndcr_ud
                b_num1 = g.dc.dcr_dc{jj,3};
                b_num2 = g.dc.dcr_dc{jj,4};
                conv_num = g.dc.dcr_dc{jj,2};
                g.dc.angdcr(jj,k) = (g.sys.theta(g.sys.bus_int(b_num1),k)-g.sys.theta(g.sys.bus_int(b_num2),k));
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
                g.dc.angdci(jj,k)=g.sys.theta(g.sys.bus_int(b_num1),k)-g.sys.theta(g.sys.bus_int(b_num2),k);
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
        
        dc_cont(0,k,10*(k-1)+1,bus_sim,flag);
        
        %% network interface for control models
        dpwf(0,k,flag);
        
        mexc_sig(k);
        smpexc(0,k,flag);
        smppi(0,k,flag);
        exc_st3(0,k,flag);
        exc_dc12(0,k,flag);
        
        mtg_sig(k);
        tg(0,k,flag);
        tg_hydro(0,k,bus_sim,flag);
        
        if g.svc.n_dcud~=0
            %% set the new line currents
            % SVC damping control...
            % Static Var Compensator
            for jj=1:g.svc.n_dcud
                l_num = g.svc.svc_dc{jj,3};
                svc_num = g.svc.svc_dc{jj,2};
                from_bus = g.sys.bus_int(line_sim(l_num,1)); 
                to_bus = g.sys.bus_int(line_sim(l_num,2));
                svc_bn = g.sys.bus_int(g.svc.svc_con(svc_num,2));
                V1 = g.sys.bus_v(from_bus,k);
                V2 = g.sys.bus_v(to_bus,k);
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
                g.tcsc.td_sig(jj,k)=abs(g.sys.bus_v(g.sys.bus_int(b_num),k));
            end
        end
        
        %% step 3b: compute dynamics and integrate
        flag = 2;
        %g.sys.sys_freq(k) = 1.0; % why?... 5/21/20 a
        % initialized as all ones on line 764
        
        mpm_sig(k);
        
        mac_ind(0,k,bus_sim,flag);
        mac_igen(0,k,bus_sim,flag);
        
        mac_sub(0,k,bus_sim,flag);
        mac_tra(0,k,bus_sim,flag);
        mac_em(0,k,bus_sim,flag);
        
        dpwf(0,k,flag);
        pss(0,k,flag);
        
        mexc_sig(k);
        smpexc(0,k,flag);
        smppi(0,k,flag);
        exc_st3(0,k,flag);
        exc_dc12(0,k,flag);
        
        mtg_sig(k);
        tg(0,k,flag);
        tg_hydro(0,k,bus_sim,flag);
        
        if g.svc.n_svc~=0
            v_svc = abs(g.sys.bus_v(g.sys.bus_int(g.svc.svc_con(:,2)),k));
            if g.svc.n_dcud~=0
                tot_states=0;
                for jj = 1:g.svc.n_dcud
                    ysvcmx = g.svc.svc_dc{jj,4};
                    ysvcmn = g.svc.svc_dc{jj,5};
                    svc_num = g.svc.svc_dc{jj,2};
                    st_state = tot_states+1; 
                    svc_states = g.svc.svc_dc{jj,6}; 
                    tot_states = tot_states+svc_states;
                    [g.svc.svc_dsig(svc_num,k),g.svc.xsvc_dc(st_state:tot_states,k),g.svc.dxsvc_dc(st_state:tot_states,k)] =...
                        svc_sud(jj,k,flag,g.svc.svc_dc{jj,1},d_sig(jj,k),ysvcmx,ysvcmn,g.svc.xsvc_dc(st_state:tot_states,k));
                end
            end
            msvc_sig(k);
            svc(0,k,bus_sim,flag,v_svc);
        end
        if g.tcsc.n_tcsc~=0
            if g.tcsc.n_tcscud~=0
                tot_states=0;
                for jj = 1:g.tcsc.n_tcscud
                    ytcscmx = g.tcsc.tcsc_dc{jj,4};
                    ytcscmn = g.tcsc.tcsc_dc{jj,5};
                    tcsc_num = g.tcsc.tcsc_dc{jj,2};
                    st_state = tot_states+1; 
                    tcsc_states = g.tcsc.tcsc_dc{jj,6}; 
                    tot_states = tot_states+tcsc_states;
                    [g.tcsc.tcsc_dsig(tcsc_num,k),g.tcsc.xtcsc_dc(st_state:tot_states,k),g.tcsc.dxtcsc_dc(st_state:tot_states,k)] =...
                        tcsc_sud(jj,k,flag,g.tcsc.tcsc_dc{jj,1},g.tcsc.td_sig(jj,k),ytcscmx,ytcscmn,g.tcsc.xtcsc_dc(st_state:tot_states,k));
                end
            end
            mtcsc_sig(k);
            tcsc(0,k,flag);
        end
        
        if g.lmod.n_lmod~=0
            ml_sig(k); 
            lmod(0,k,flag); 
        end
        
        if g.rlmod.n_rlmod~=0
            rml_sig(k);
            rlmod(0,k,flag);
        end
        
        %% pwrmod - copied from v2.3 - 06/01/20 -thad
        if g.pwr.n_pwrmod~=0
            Pst = cell(g.pwr.n_pwrmod,1);
            Qst = Pst;
            for index=1:g.pwr.n_pwrmod
                Pst{index} = pwrmod_p_sigst{index}(:,k);
                Qst{index} = pwrmod_q_sigst{index}(:,k);
            end
            [~,~,dp,dq,~,~] = pwrmod_dyn(Pst,Qst,bus,t,k,flag,g.pwr.n_pwrmod);
            if (~iscell(dp) || ~iscell(dq)) 
                error('Error in pwrmod_dyn, dp and dq must be cells'); 
            end
            if (any(size(dp)-[g.pwr.n_pwrmod 1]) || any(size(dq)-[g.pwr.n_pwrmod 1]))
                error('Dimension error in pwrmod_dyn'); 
            end
            for index=1:g.pwr.n_pwrmod
                if ((size(dp{index},2)~=1) || (size(dq{index},2)~=1))
                    error('Dimension error in pwrmod_dyn'); 
                end
                if size(dp{index},1)~=size(dpwrmod_p_sigst{index},1)
                    error('Dimension error in pwrmod_dyn'); 
                end
                if size(dq{index},1)~=size(dpwrmod_q_sigst{index},1)
                    error('Dimension error in pwrmod_dyn'); 
                end
                dpwrmod_p_sigst{index}(:,k) = dp{index};
                dpwrmod_q_sigst{index}(:,k) = dq{index};
            end
            [P,Q,~,~] = pwrmod_dyn(Pst,Qst,bus,t,k,1,g.pwr.n_pwrmod); %update pwrmod_p_sig and pwrmod_q_sig
            if (length(P)~=g.pwr.n_pwrmod) || (length(Q)~=g.pwr.n_pwrmod)
                error('Dimension error in pwrmod_dyn'); 
            end
            g.pwr.pwrmod_p_sig(:,k) = P;
            g.pwr.pwrmod_q_sig(:,k) = Q;
            pwrmod_p(0,k,bus_sim,flag);
            pwrmod_q(0,k,bus_sim,flag);
            clear P Q Pst Qst dp dq index
        end
        
        %% ivm modulation - copied from v2.3 - 06/01/20 -thad
        if n_ivm>0
            dst = cell(n_ivm,1);
            est = dst;
            for index=1:n_ivm
                dst{index} = ivmmod_d_sigst{index}(:,k);
                est{index} = ivmmod_e_sigst{index}(:,k);
            end
            [d,e,~,~,~,~] = ivmmod_dyn(dst,est,bus,t,k,1); %get internal voltage signals
            if (length(d)~=n_ivm) || (length(e)~=n_ivm)
                error('Dimension error in ivmmod_dyn');
            end
            ivmmod_d_sig(:,k) = d;
            ivmmod_e_sig(:,k) = e;
            mac_ivm(0,k,bus_sim,flag);
            [~,~,dd,de,~,~] = ivmmod_dyn(dst,est,bus,t,k,flag);
            if (~iscell(dd) || ~iscell(de))
                error('Error in ivmmod_dyn, dd and de must be cells'); 
            end
            if (any(size(dd)-[n_ivm 1]) || any(size(de)-[n_ivm 1]))
                error('Dimension error in ivmmod_dyn'); 
            end
            for index=1:n_ivm
                if ((size(dd{index},2)~=1) || (size(de{index},2)~=1))
                    error('Dimension error in ivmmod_dyn'); 
                end
                if size(dd{index},1)~=size(divmmod_d_sigst{index},1)
                    error('Dimension error in ivmmod_dyn'); 
                end
                if size(de{index},1)~=size(divmmod_e_sigst{index},1)
                    error('Dimension error in ivmmod_dyn'); 
                end
                divmmod_d_sigst{index}(:,k) = dd{index};
                divmmod_e_sigst{index}(:,k) = de{index};
            end
            clear d e dd de dst est
        end
        
        
        %% integrate dc at ten times rate (DC Stuff? 5/14/20)
        mdc_sig(k);
        if g.dc.n_conv~=0
            hdc_sol = h_sol/10;
            for kk = 1:10
                kdc=10*(k-1)+kk;
                [g.dc.xdcr_dc(:,kdc:kdc+1),g.dc.dxdcr_dc(:,kdc:kdc+1),g.dc.xdci_dc(:,kdc:kdc+1),g.dc.dxdci_dc(:,kdc:kdc+1)] = ...
                    dc_sim(k,kk,g.dc.dcr_dc,g.dc.dci_dc,g.dc.xdcr_dc(:,kdc),g.dc.xdci_dc(:,kdc),bus_sim,hdc_sol); % dc_sim 
            end
        else
            dc_cont(0,k,k,bus_sim,2);
            dc_line(0,k,k,bus_sim,2);
        end
        
        %% following statements are predictor steps
        j = k+1;
        g.mac.mac_ang(:,j) = g.mac.mac_ang(:,k) + h_sol*g.mac.dmac_ang(:,k);
        g.mac.mac_spd(:,j) = g.mac.mac_spd(:,k) + h_sol*g.mac.dmac_spd(:,k);
        g.mac.edprime(:,j) = g.mac.edprime(:,k) + h_sol*g.mac.dedprime(:,k);
        g.mac.eqprime(:,j) = g.mac.eqprime(:,k) + h_sol*g.mac.deqprime(:,k);
        g.mac.psikd(:,j) = g.mac.psikd(:,k) + h_sol*g.mac.dpsikd(:,k);
        g.mac.psikq(:,j) = g.mac.psikq(:,k) + h_sol*g.mac.dpsikq(:,k);
        
        % Exciter integration
        g.exc.Efd(:,j) = g.exc.Efd(:,k) + h_sol*g.exc.dEfd(:,k);
        g.exc.V_R(:,j) = g.exc.V_R(:,k) + h_sol*g.exc.dV_R(:,k);
        g.exc.V_As(:,j) = g.exc.V_As(:,k) + h_sol*g.exc.dV_As(:,k);
        g.exc.R_f(:,j) = g.exc.R_f(:,k) + h_sol*g.exc.dR_f(:,k);
        g.exc.V_TR(:,j) = g.exc.V_TR(:,k) + h_sol*g.exc.dV_TR(:,k);
        
        if n_dpw ~= 0
            % only calculate if dpw filter is used
            sdpw1(:,j) = sdpw1(:,k) + h_sol*dsdpw1(:,k);
            sdpw2(:,j) = sdpw2(:,k) + h_sol*dsdpw2(:,k);
            sdpw3(:,j) = sdpw3(:,k) + h_sol*dsdpw3(:,k);
            sdpw4(:,j) = sdpw4(:,k) + h_sol*dsdpw4(:,k);
            sdpw5(:,j) = sdpw5(:,k) + h_sol*dsdpw5(:,k);
            sdpw6(:,j) = sdpw6(:,k) + h_sol*dsdpw6(:,k);
        end
        
        g.pss.pss1(:,j) = g.pss.pss1(:,k) + h_sol*g.pss.dpss1(:,k);
        g.pss.pss2(:,j) = g.pss.pss2(:,k) + h_sol*g.pss.dpss2(:,k);
        g.pss.pss3(:,j) = g.pss.pss3(:,k) + h_sol*g.pss.dpss3(:,k);
        
        % modified to g - thad
        g.tg.tg1(:,j) = g.tg.tg1(:,k) + h_sol*g.tg.dtg1(:,k);
        g.tg.tg2(:,j) = g.tg.tg2(:,k) + h_sol*g.tg.dtg2(:,k);
        g.tg.tg3(:,j) = g.tg.tg3(:,k) + h_sol*g.tg.dtg3(:,k);
        g.tg.tg4(:,j) = g.tg.tg4(:,k) + h_sol*g.tg.dtg4(:,k);
        g.tg.tg5(:,j) = g.tg.tg5(:,k) + h_sol*g.tg.dtg5(:,k);
        
        % induction motor integrations
        if g.ind.n_mot ~= 0
            g.ind.vdp(:,j) = g.ind.vdp(:,k) + h_sol*g.ind.dvdp(:,k);
            g.ind.vqp(:,j) = g.ind.vqp(:,k) + h_sol*g.ind.dvqp(:,k);
            g.ind.slip(:,j) = g.ind.slip(:,k) + h_sol*g.ind.dslip(:,k);
        end
        
        % induction generator integrations
        if g.igen.n_ig ~=0
            g.igen.vdpig(:,j) = g.igen.vdpig(:,k) + h_sol*g.igen.dvdpig(:,k);
            g.igen.vqpig(:,j) = g.igen.vqpig(:,k) + h_sol*g.igen.dvqpig(:,k);
            g.igen.slig(:,j) = g.igen.slig(:,k) + h_sol*g.igen.dslig(:,k);
        end
        
        % svc
        if g.svc.n_svc ~= 0
            g.svc.B_cv(:,j) = g.svc.B_cv(:,k) + h_sol*g.svc.dB_cv(:,k);
            g.svc.B_con(:,j) = g.svc.B_con(:,k) + h_sol*g.svc.dB_con(:,k);
            g.svc.xsvc_dc(:,j) = g.svc.xsvc_dc(:,k) + h_sol* g.svc.dxsvc_dc(:,k);
        end
        
        %tcsc
        if g.tcsc.n_tcsc ~= 0
            g.tcsc.B_tcsc(:,j) = g.tcsc.B_tcsc(:,k) + h_sol*g.tcsc.dB_tcsc(:,k);
            g.tcsc.xtcsc_dc(:,j) = g.tcsc.xtcsc_dc(:,k) + h_sol* g.tcsc.dxtcsc_dc(:,k);
        end
        
        if g.lmod.n_lmod~=0
            g.lmod.lmod_st(:,j) = g.lmod.lmod_st(:,k) + h_sol*g.lmod.dlmod_st(:,k); % line using g
        end
        
        if g.rlmod.n_rlmod~=0
            g.rlmod.rlmod_st(:,j) = g.rlmod.rlmod_st(:,k)+h_sol*g.rlmod.drlmod_st(:,k);
        end
        %% Copied from v2.3 - 06/01/20 - thad
        g.pwr.pwrmod_p_st(:,j) = g.pwr.pwrmod_p_st(:,k)+h_sol*g.pwr.dpwrmod_p_st(:,k);
        g.pwr.pwrmod_q_st(:,j) = g.pwr.pwrmod_q_st(:,k)+h_sol*g.pwr.dpwrmod_q_st(:,k);
        %% pwrmod
        if g.pwr.n_pwrmod~=0
            for index=1:g.pwr.n_pwrmod
                pwrmod_p_sigst{index}(:,j) = pwrmod_p_sigst{index}(:,k)+h_sol*dpwrmod_p_sigst{index}(:,k);
                pwrmod_q_sigst{index}(:,j) = pwrmod_q_sigst{index}(:,k)+h_sol*dpwrmod_q_sigst{index}(:,k);
            end
        end
        %% ivmmod
        if n_ivm~=0
            for index=1:n_ivm
                ivmmod_d_sigst{index}(:,j) = ivmmod_d_sigst{index}(:,k)+h_sol*divmmod_d_sigst{index}(:,k);
                ivmmod_e_sigst{index}(:,j) = ivmmod_e_sigst{index}(:,k)+h_sol*divmmod_e_sigst{index}(:,k);
            end
        end
        
        
        %% Flag = 1
        % begining of solutions as j 
        flag = 1;
        % mach_ref(j) = mac_ang(syn_ref,j);
        g.sys.mach_ref(j) = 0;
        % perform network interface calculations again with predicted states
        mpm_sig(j);
        mac_ind(0,j,bus_sim,flag);
        mac_igen(0,j,bus_sim,flag);
        mac_sub(0,j,bus_sim,flag);
        mac_tra(0,j,bus_sim,flag);
        mac_em(0,j,bus_sim,flag);
        mac_ivm(0,j,bus_sim,flag); 
        
        % assume Vdc remains unchanged for first pass through dc controls interface
        mdc_sig(j);
        dc_cont(0,j,10*(j-1)+1,bus_sim,flag);
        
        % Calculate current injections and bus voltages and angles
        if j >= sum(k_inc(1:3))+1
            % fault cleared
            bus_sim = bus_pf2;
            g.sys.bus_int = bus_intpf2;
            Y1 = Y_gpf2;
            Y2 = Y_gncpf2;
            Y3 = Y_ncgpf2;
            Y4 = Y_ncpf2;
            Vr1 = V_rgpf2;
            Vr2 = V_rncpf2;
            bo = bopf2;
            %h_sol = i_simu(j,ks,k_inc,h,bus_sim,Y1,Y2,Y3,Y4,Vr1,Vr2,bo);
        elseif j >=sum(k_inc(1:2))+1
            % near bus cleared
            bus_sim = bus_pf1;
            g.sys.bus_int = bus_intpf1;
            Y1 = Y_gpf1;
            Y2 = Y_gncpf1;
            Y3 = Y_ncgpf1;
            Y4 = Y_ncpf1;
            Vr1 = V_rgpf1;
            Vr2 = V_rncpf1;
            bo = bopf1;
            %h_sol = i_simu(j,ks,k_inc,h,bus_sim,Y1,Y2,Y3,Y4,Vr1,Vr2,bo);
        elseif j>=k_inc(1)+1
            % fault applied
            bus_sim = bus_f;
            g.sys.bus_int = bus_intf;
            Y1 = Y_gf;
            Y2 = Y_gncf;
            Y3 = Y_ncgf;
            Y4 = Y_ncf;
            Vr1 = V_rgf;
            Vr2 = V_rncf;
            bo = bof;
            %h_sol = i_simu(j,ks,k_inc,h,bus_sim,Y1,Y2,Y3,Y4,Vr1,Vr2,bo);
        elseif k<k_inc(1)+1  % JHC - DKF thinks k should be j
            % pre fault
            bus_sim = bus;
            g.sys.bus_int = bus_intprf;
            Y1 = Y_gprf;
            Y2 = Y_gncprf;
            Y3 = Y_ncgprf;
            Y4 = Y_ncprf;
            Vr1 = V_rgprf;
            Vr2 = V_rncprf;
            bo = boprf;
            %h_sol = i_simu(j,ks,k_inc,h,bus_sim,Y1,Y2,Y3,Y4,Vr1,Vr2,bo);
        end
        
        % apply gen trip - copied from v2.3 - 06/01/20 - thad
        if sum(mac_trip_flags)>0.5
            genBuses = g.mac.mac_con(mac_trip_flags==1,2);
            for kB=1:length(genBuses)
                nL = find(genBuses(kB)==line_sim(:,1) | genBuses(kB)==line_sim(:,2));
                if isempty(nL)
                    error('nL is empty.'); 
                end
                line_sim(nL,4) = 1e7; %make reactance infinity
            end
            [Y1,Y2,Y3,Y4,Vr1,Vr2,bo] = red_ybus(bus_sim,line_sim);
            clear nL kB genBuses
        end
        
        %% solve
        h_sol = i_simu(j,ks,k_inc,h,bus_sim,Y1,Y2,Y3,Y4,Vr1,Vr2,bo);
        
        g.mac.vex(:,j) = g.mac.vex(:,k);
        g.dc.cur_ord(:,j) = g.dc.cur_ord(:,k);
        % calculate the new value of bus angles rectifier user defined control
        if g.dc.ndcr_ud~=0
            tot_states=0;
            for jj = 1:g.dc.ndcr_ud
                b_num1 = g.dc.dcr_dc{jj,3};
                b_num2 = g.dc.dcr_dc{jj,4};
                conv_num = g.dc.dcr_dc{jj,2};
                g.dc.angdcr(jj,j) = g.sys.theta(g.sys.bus_int(b_num1),j)-g.sys.theta(g.sys.bus_int(b_num2),j);
                g.dc.dcrd_sig(jj,j)=g.dc.angdcr(jj,j);
                st_state = tot_states+1; 
                dcr_states = g.dc.dcr_dc{jj,7}; 
                tot_states = tot_states+dcr_states;
                ydcrmx = g.dc.dcr_dc{jj,5};
                ydcrmn = g.dc.dcr_dc{jj,6};
                g.dc.dcr_dsig(jj,j) = ...
                    dcr_sud(jj,j,flag,g.dc.dcr_dc{jj,1},g.dc.dcrd_sig(jj,j),ydcrmx,ydcrmn,g.dc.xdcr_dc(st_state:tot_states,10*(j-1)+1));
            end
        end
        if g.dc.ndci_ud~=0
            % calculate the new value of bus angles inverter user defined control
            for jj = 1:g.dc.ndci_ud
                tot_states=0;
                b_num1 = g.dc.dci_dc{jj,3};
                b_num2 = g.dc.dci_dc{jj,4};
                conv_num = g.dc.dci_dc{jj,2};
                g.dc.angdci(jj,j) = g.sys.theta(g.sys.bus_int(b_num1),j)-g.sys.theta(g.sys.bus_int(b_num2),j);
                g.dc.dcid_sig(jj,j) = (g.dc.angdci(jj,j)-g.dc.angdci(jj,k))/(t(j)-t(k));
                st_state = tot_states+1; 
                dci_states = g.dc.dci_dc{jj,7}; 
                tot_states = tot_states+dci_states;
                ydcimx = g.dc.dci_dc{jj,5};
                ydcimn = g.dc.dci_dc{jj,6};
                g.dc.dci_dsig(jj,j) = ...
                    dci_sud(jj,j,flag,g.dc.dci_dc{jj,1},g.dc.dcid_sig(jj,j),ydcimx,ydcimn,g.dc.xdci_dc(st_state:tot_states,10*(j-1)+1));
            end
        end
        
        %% network interface for control models - 'corrector' step
        dc_cont(0,j,10*(j-1)+1,bus_sim,flag);
        
        dpwf(0,j,flag);
        pss(0,j,flag);
        
        mexc_sig(j); % modulation
        smpexc(0,j,flag);
        smppi(0,j,flag);
        exc_st3(0,j,flag);
        exc_dc12(0,j,flag);
        
        tg(0,j,flag);
        tg_hydro(0,j,bus_sim,flag);
        
        if g.svc.n_dcud~=0
            % set the new line currents
            for jj=1:g.svc.n_dcud
                l_num = g.svc.svc_dc{jj,3};svc_num = g.svc.svc_dc{jj,2};
                from_bus = g.sys.bus_int(line_sim(l_num,1)); 
                to_bus = g.sys.bus_int(line_sim(l_num,2));
                svc_bn = g.sys.bus_int(g.svc.svc_con(svc_num,2));
                V1 = g.sys.bus_v(from_bus,j);
                V2 = g.sys.bus_v(to_bus,j);
                R = line_sim(l_num,3);
                X = line_sim(l_num,4);
                B = line_sim(l_num,5);
                g.dc.tap = line_sim(l_num,6);phi = line_sim(l_num,7);
                [l_if,l_it] = line_cur(V1,V2,R,X,B,g.dc.tap,phi);
                if svc_bn == from_bus
                    d_sig(jj,j)=abs(l_if);
                elseif svc_bn==to_bus
                    d_sig(jj,j)=abs(l_it);
                end
            end
        end
        
        if g.tcsc.n_tcscud~=0
            % set the new line currents
            for jj=1:g.tcsc.n_tcscud
                b_num = g.tcsc.tcsc_dc{jj,3};
                tcsc_num = g.tcsc.tcsc_dc{jj,2};
                g.tcsc.td_sig(jj,j) = abs(g.sys.bus_v(g.sys.bus_int(b_num),j));
            end
        end
        
        %% Flag = 2, for 'corrector step' d's
        flag = 2;
        mac_ind(0,j,bus_sim,flag);
        mac_igen(0,j,bus_sim,flag);
        mac_sub(0,j,bus_sim,flag);
        mac_tra(0,j,bus_sim,flag);
        mac_em(0,j,bus_sim,flag);
        
        dpwf(0,j,flag);
        pss(0,j,flag);
        
        mexc_sig(j); % modulation
        smpexc(0,j,flag);
        smppi(0,j,flag);
        exc_st3(0,j,flag);
        exc_dc12(0,j,flag);
        
        mtg_sig(j);% modulation
        tg(0,j,flag);
        tg_hydro(0,j,bus_sim,flag);
        
        if g.svc.n_svc~=0
            msvc_sig(j);% modulation
            if g.svc.n_dcud~=0
                tot_states=0;
                for jj = 1:g.svc.n_dcud
                    ysvcmx = g.svc.svc_dc{jj,4};
                    ysvcmn = g.svc.svc_dc{jj,5};
                    svc_num = g.svc.svc_dc{jj,2};
                    st_state = tot_states+1; 
                    svc_states = g.svc.svc_dc{jj,6}; 
                    tot_states = tot_states+svc_states;
                    [g.svc.svc_dsig(svc_num,j),g.svc.xsvc_dc(st_state:tot_states,j),g.svc.dxsvc_dc(st_state:tot_states,j)] =...
                        svc_sud(jj,j,flag,g.svc.svc_dc{jj,1},d_sig(jj,j),ysvcmx,ysvcmn,g.svc.xsvc_dc(st_state:tot_states,j));
                end
            end
            v_svc = abs(g.sys.bus_v(g.sys.bus_int(g.svc.svc_con(:,2)),j));
            bus_sim = svc(0,j,bus_sim,flag,v_svc);
        end
        
        if g.tcsc.n_tcsc~=0
            mtcsc_sig(j); % this has changed since v 2.3... % modulation
            if g.tcsc.n_tcscud~=0
                tot_states=0;
                for jj = 1:g.tcsc.n_tcscud
                    ytcscmx = g.tcsc.tcsc_dc{jj,4};
                    ytcscmn = g.tcsc.tcsc_dc{jj,5};
                    tcsc_num = g.tcsc.tcsc_dc{jj,2};
                    st_state = tot_states+1; 
                    tcsc_states = g.tcsc.tcsc_dc{jj,6}; 
                    tot_states = tot_states+tcsc_states;
                    [g.tcsc.tcsc_dsig(tcsc_num,j),g.tcsc.xtcsc_dc(st_state:tot_states,j),g.tcsc.dxtcsc_dc(st_state:tot_states,j)] =...
                        tcsc_sud(jj,j,flag,g.tcsc.tcsc_dc{jj,1},g.tcsc.td_sig(jj,j),ytcscmx,ytcscmn,g.tcsc.xtcsc_dc(st_state:tot_states,j));
                end
            end
            tcsc(0,j,flag);
        end
        
        % modified to handle g - thad 06/01/20
        if g.lmod.n_lmod~=0
            ml_sig(j); % removed t - thad % modulation
            lmod(0,j,flag); % removed bus - thad
        end
        if g.rlmod.n_rlmod~=0
            rml_sig(j); % modulation
            rlmod(0,j,flag);
        end
        
        % copied from v2.3 - thad - 06/01/20
        if g.pwr.n_pwrmod~=0
            Pst = cell(g.pwr.n_pwrmod,1);
            Qst = Pst;
            for index=1:g.pwr.n_pwrmod
                Pst{index} = pwrmod_p_sigst{index}(:,j);
                Qst{index} = pwrmod_q_sigst{index}(:,j);
            end
            [~,~,dp,dq,~,~] = pwrmod_dyn(Pst,Qst,bus,t,j,flag,g.pwr.n_pwrmod);
            if (~iscell(dp) || ~iscell(dq))
                error('Error in pwrmod_dyn, dp and dq must be cells'); 
            end
            if (any(size(dp)-[g.pwr.n_pwrmod 1]) || any(size(dq)-[g.pwr.n_pwrmod 1]))
                error('Dimension error in pwrmod_dyn'); 
            end
            
            for index=1:g.pwr.n_pwrmod
                if ((size(dp{index},2)~=1) || (size(dq{index},2)~=1))
                    error('Dimension error in pwrmod_dyn'); 
                end
                if size(dp{index},1)~=size(dpwrmod_p_sigst{index},1)
                    error('Dimension error in pwrmod_dyn'); 
                end
                if size(dq{index},1)~=size(dpwrmod_q_sigst{index},1)
                    error('Dimension error in pwrmod_dyn'); 
                end
                dpwrmod_p_sigst{index}(:,j) = dp{index};
                dpwrmod_q_sigst{index}(:,j) = dq{index};
            end
            [P,Q,~,~,~,~] = pwrmod_dyn(Pst,Qst,bus,t,j,1,g.pwr.n_pwrmod); %update pwrmod_p_sig and pwrmod_q_sig
            if (length(P)~=g.pwr.n_pwrmod) || (length(Q)~=g.pwr.n_pwrmod)
                error('Dimension error in pwrmod_dyn');
            end
            g.pwr.pwrmod_p_sig(:,j) = P;
            g.pwr.pwrmod_q_sig(:,j) = Q;
            pwrmod_p(0,j,bus_sim,flag);
            pwrmod_q(0,j,bus_sim,flag);
            clear P Q Pst Qst dp dq index
        end
        
        if n_ivm>0
            dst = cell(n_ivm,1);
            est = dst;
            for index=1:n_ivm
                dst{index} = ivmmod_d_sigst{index}(:,j);
                est{index} = ivmmod_e_sigst{index}(:,j);
            end
            [d,e,~,~,~,~] = ivmmod_dyn(dst,est,bus,t,j,1);
            if (length(d)~=n_ivm) || (length(e)~=n_ivm)
                error('Dimension error in ivmmod_dyn'); 
            end
            ivmmod_d_sig(:,j) = d;
            ivmmod_e_sig(:,j) = e;
            mac_ivm(0,j,bus_sim,flag);
            [~,~,dd,de,~,~] = ivmmod_dyn(dst,est,bus,t,j,flag);
            if (~iscell(dd) || ~iscell(de))
                error('Error in ivmmod_dyn, dd and de must be cells'); 
            end
            if (any(size(dd)-[n_ivm 1]) || any(size(de)-[n_ivm 1]))
                error('Dimension error in ivmmod_dyn'); 
            end
            for index=1:n_ivm
                if ((size(dd{index},2)~=1) || (size(de{index},2)~=1))
                    error('Dimension error in ivmmod_dyn'); 
                end
                if size(dd{index},1)~=size(divmmod_d_sigst{index},1)
                    error('Dimension error in ivmmod_dyn'); 
                end
                if size(de{index},1)~=size(divmmod_e_sigst{index},1)
                    error('Dimension error in ivmmod_dyn'); 
                end
                divmmod_d_sigst{index}(:,j) = dd{index};
                divmmod_e_sigst{index}(:,j) = de{index};
            end
            clear d e dd de dst est
        end
        % end copied from...
        
        %%integrate dc at ten times rate (DC Stuff? 6/09/20)
        if g.dc.n_conv~=0
            hdc_sol = h_sol/10;
            for kk = 1:10
                jdc=10*(j-1)+kk;
                [g.dc.xdcr_dc(:,jdc:jdc+1),g.dc.dxdcr_dc(:,jdc:jdc+1),g.dc.xdci_dc(:,jdc:jdc+1),g.dc.dxdci_dc(:,jdc:jdc+1)] = ...
                    dc_sim(j,kk,g.dc.dcr_dc,g.dc.dci_dc,g.dc.xdcr_dc(:,jdc),g.dc.xdci_dc(:,jdc),bus_sim,hdc_sol);
            end
        else
            dc_cont(0,j,j,bus_sim,2);
            dc_line(0,j,j,bus_sim,2);
        end
        
        %% following statements are corrector steps (Actual RK2 computation)
        g.mac.mac_ang(:,j) = g.mac.mac_ang(:,k) + h_sol*(g.mac.dmac_ang(:,k)+g.mac.dmac_ang(:,j))/2.;
        g.mac.mac_spd(:,j) = g.mac.mac_spd(:,k) + h_sol*(g.mac.dmac_spd(:,k)+g.mac.dmac_spd(:,j))/2.;
        g.mac.edprime(:,j) = g.mac.edprime(:,k) + h_sol*(g.mac.dedprime(:,k)+g.mac.dedprime(:,j))/2.;
        g.mac.eqprime(:,j) = g.mac.eqprime(:,k) + h_sol*(g.mac.deqprime(:,k)+g.mac.deqprime(:,j))/2.;
        g.mac.psikd(:,j) = g.mac.psikd(:,k) + h_sol*(g.mac.dpsikd(:,k)+g.mac.dpsikd(:,j))/2.;
        g.mac.psikq(:,j) = g.mac.psikq(:,k) + h_sol*(g.mac.dpsikq(:,k)+g.mac.dpsikq(:,j))/2.;
        
        % exciter integration
        g.exc.Efd(:,j) = g.exc.Efd(:,k) + h_sol*(g.exc.dEfd(:,k)+g.exc.dEfd(:,j))/2.;
        g.exc.V_R(:,j) = g.exc.V_R(:,k) + h_sol*(g.exc.dV_R(:,k)+g.exc.dV_R(:,j))/2.;
        g.exc.V_As(:,j) = g.exc.V_As(:,k) + h_sol*(g.exc.dV_As(:,k)+g.exc.dV_As(:,j))/2.;
        g.exc.R_f(:,j) = g.exc.R_f(:,k) + h_sol*(g.exc.dR_f(:,k)+g.exc.dR_f(:,j))/2.;
        g.exc.V_TR(:,j) = g.exc.V_TR(:,k) + h_sol*(g.exc.dV_TR(:,k)+g.exc.dV_TR(:,j))/2.;
        
        % removed extra 1 in global names. - thad 07/06/20
        if n_dpw ~= 0
            % only calculate if dpw filter is used
            sdpw1(:,j) = sdpw1(:,k) +h_sol*(dsdpw1(:,k)+dsdpw1(:,j))/2.;
            sdpw2(:,j) = sdpw2(:,k) +h_sol*(dsdpw2(:,k)+dsdpw2(:,j))/2.;
            sdpw3(:,j) = sdpw3(:,k) +h_sol*(dsdpw3(:,k)+dsdpw3(:,j))/2.;
            sdpw4(:,j) = sdpw4(:,k) +h_sol*(dsdpw4(:,k)+dsdpw4(:,j))/2.;
            sdpw5(:,j) = sdpw5(:,k) +h_sol*(dsdpw5(:,k)+dsdpw5(:,j))/2.;
            sdpw6(:,j) = sdpw6(:,k) +h_sol*(dsdpw6(:,k)+dsdpw6(:,j))/2.;
        end
        
        g.pss.pss1(:,j) = g.pss.pss1(:,k) +h_sol*(g.pss.dpss1(:,k)+g.pss.dpss1(:,j))/2.;
        g.pss.pss2(:,j) = g.pss.pss2(:,k) +h_sol*(g.pss.dpss2(:,k)+g.pss.dpss2(:,j))/2.;
        g.pss.pss3(:,j) = g.pss.pss3(:,k) +h_sol*(g.pss.dpss3(:,k)+g.pss.dpss3(:,j))/2.;
        
        % modified to g
        g.tg.tg1(:,j) = g.tg.tg1(:,k) + h_sol*(g.tg.dtg1(:,k) + g.tg.dtg1(:,j))/2.;
        g.tg.tg2(:,j) = g.tg.tg2(:,k) + h_sol*(g.tg.dtg2(:,k) + g.tg.dtg2(:,j))/2.;
        g.tg.tg3(:,j) = g.tg.tg3(:,k) + h_sol*(g.tg.dtg3(:,k) + g.tg.dtg3(:,j))/2.;
        g.tg.tg4(:,j) = g.tg.tg4(:,k) + h_sol*(g.tg.dtg4(:,k) + g.tg.dtg4(:,j))/2.;
        g.tg.tg5(:,j) = g.tg.tg5(:,k) + h_sol*(g.tg.dtg5(:,k) + g.tg.dtg5(:,j))/2.;
        
        % induction motor integrations
        if g.ind.n_mot ~= 0
            g.ind.vdp(:,j) = g.ind.vdp(:,k) + h_sol*(g.ind.dvdp(:,j) + g.ind.dvdp(:,k))/2.;
            g.ind.vqp(:,j) = g.ind.vqp(:,k) + h_sol*(g.ind.dvqp(:,j) + g.ind.dvqp(:,k))/2.;
            g.ind.slip(:,j) = g.ind.slip(:,k) + h_sol*(g.ind.dslip(:,j) + g.ind.dslip(:,k))/2.;
        end
        
        % induction generator integrations
        if g.igen.n_ig ~=0
            g.igen.vdpig(:,j) = g.igen.vdpig(:,k) + h_sol*(g.igen.dvdpig(:,j) + g.igen.dvdpig(:,k))/2.;
            g.igen.vqpig(:,j) = g.igen.vqpig(:,k) + h_sol*(g.igen.dvqpig(:,j) + g.igen.dvqpig(:,k))/2.;
            g.igen.slig(:,j) = g.igen.slig(:,k) + h_sol*(g.igen.dslig(:,j) + g.igen.dslig(:,k))/2.;
        end
        
        % svc
        if g.svc.n_svc ~= 0
            g.svc.B_cv(:,j) = g.svc.B_cv(:,k) + h_sol*(g.svc.dB_cv(:,j) + g.svc.dB_cv(:,k))/2.;
            g.svc.B_con(:,j) = g.svc.B_con(:,k) + h_sol*(g.svc.dB_con(:,j) + g.svc.dB_con(:,k))/2.;
            g.svc.xsvc_dc(:,j) = g.svc.xsvc_dc(:,k) + h_sol*(g.svc.dxsvc_dc(:,j) + g.svc.dxsvc_dc(:,k))/2.;
        end
        
        %tcsc
        if g.tcsc.n_tcsc ~= 0
            g.tcsc.B_tcsc(:,j) = g.tcsc.B_tcsc(:,k) + h_sol*(g.tcsc.dB_tcsc(:,j) + g.tcsc.dB_tcsc(:,k))/2.;
            g.tcsc.xtcsc_dc(:,j) = g.tcsc.xtcsc_dc(:,k) + h_sol*(g.tcsc.dxtcsc_dc(:,j) + g.tcsc.dxtcsc_dc(:,k))/2.;  
        end
        
        if g.lmod.n_lmod~=0
            g.lmod.lmod_st(:,j) = g.lmod.lmod_st(:,k) + h_sol*(g.lmod.dlmod_st(:,j) + g.lmod.dlmod_st(:,k))/2.; % modified line with g
        end
        if g.rlmod.n_rlmod~=0
            g.rlmod.rlmod_st(:,j) = g.rlmod.rlmod_st(:,k) + h_sol*(g.rlmod.drlmod_st(:,j) + g.rlmod.drlmod_st(:,k))/2.;
        end
        
        % Copied from v2.3 - 06/01/20 - thad
        g.pwr.pwrmod_p_st(:,j) = g.pwr.pwrmod_p_st(:,k)+h_sol*(g.pwr.dpwrmod_p_st(:,j) + g.pwr.dpwrmod_p_st(:,k))/2;
        g.pwr.pwrmod_q_st(:,j) = g.pwr.pwrmod_q_st(:,k)+h_sol*(g.pwr.dpwrmod_q_st(:,j) + g.pwr.dpwrmod_q_st(:,k))/2;
        if g.pwr.n_pwrmod~=0
            for index=1:g.pwr.n_pwrmod
                % make global? -thad 07/06/20
                pwrmod_p_sigst{index}(:,j) = pwrmod_p_sigst{index}(:,k)+h_sol*(dpwrmod_p_sigst{index}(:,j) + dpwrmod_p_sigst{index}(:,k))/2;
                pwrmod_q_sigst{index}(:,j) = pwrmod_q_sigst{index}(:,k)+h_sol*(dpwrmod_q_sigst{index}(:,j) + dpwrmod_q_sigst{index}(:,k))/2;
            end
        end
        
        if n_ivm~=0
            for index=1:n_ivm
                % make global? -thad 07/06/20
                ivmmod_d_sigst{index}(:,j) = ivmmod_d_sigst{index}(:,k)+h_sol*(divmmod_d_sigst{index}(:,j) + divmmod_d_sigst{index}(:,k))/2;
                ivmmod_e_sigst{index}(:,j) = ivmmod_e_sigst{index}(:,k)+h_sol*(divmmod_e_sigst{index}(:,j) + divmmod_e_sigst{index}(:,k))/2;
            end
        end
        
        %% Live plot call
        if g.sys.livePlotFlag
           livePlot
        end
        
    end
    % counter increment
    kt = kt + k_inc(ks);
    ks = ks+1;
end% end simulation loop

% calculation of line currents post sim
V1 = g.sys.bus_v(g.sys.bus_int(line(:,1)),:);
V2 = g.sys.bus_v(g.sys.bus_int(line(:,2)),:);
R = line(:,3); 
X = line(:,4); 
B = line(:,5);
g.dc.tap = line(:,6); 
phi = line(:,7);

[ilf,ilt]=line_cur(V1,V2,R,X,B,g.dc.tap,phi);%line currents
[sInjF,sInjT]=line_pq(V1,V2,R,X,B,g.dc.tap,phi);% 'line flows' - complex power injection at bus

% full sim timing
et = toc;
ets = num2str(et);
g.sys.ElapsedNonLinearTime = ets;
disp(['elapsed time = ' ets 's'])
disp('*** End simulation.')
disp(' ')

%% Clean up logged DC variables to length of DC simulated time.
g.dc.t_dc = g.dc.t_dc(1:length(g.dc.t_dc)-10);
g.dc.i_dc = g.dc.i_dc(:,1:length(g.dc.t_dc));
g.dc.i_dcr = g.dc.i_dcr(:,1:length(g.dc.t_dc));
g.dc.i_dci = g.dc.i_dci(:,1:length(g.dc.t_dc));
g.dc.alpha = g.dc.alpha(:,1:length(g.dc.t_dc));
g.dc.gamma = g.dc.gamma(:,1:length(g.dc.t_dc));
g.dc.Vdc = g.dc.Vdc(:,1:length(g.dc.t_dc));
g.dc.v_conr = g.dc.v_conr(:,1:length(g.dc.t_dc));
g.dc.v_coni = g.dc.v_coni(:,1:length(g.dc.t_dc));
g.dc.dv_conr = g.dc.dv_conr(:,1:length(g.dc.t_dc));
g.dc.dv_coni = g.dc.dv_coni(:,1:length(g.dc.t_dc));
g.dc.xdcr_dc = g.dc.xdcr_dc(:,1:length(g.dc.t_dc));
g.dc.xdci_dc = g.dc.xdci_dc(:,1:length(g.dc.t_dc));
g.dc.dxdcr_dc = g.dc.dxdcr_dc(:,1:length(g.dc.t_dc));
g.dc.dxdci_dc = g.dc.dxdci_dc(:,1:length(g.dc.t_dc));
g.dc.dv_dcc = g.dc.dv_dcc(:,1:length(g.dc.t_dc));
g.dc.v_dcc = g.dc.v_dcc(:,1:length(g.dc.t_dc));
g.dc.di_dci = g.dc.di_dci(:,1:length(g.dc.t_dc));
g.dc.di_dcr = g.dc.di_dcr(:,1:length(g.dc.t_dc));