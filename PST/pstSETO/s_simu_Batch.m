%S_SIMU_BATCH runs a non-linear simulation of a given data file.
% S_SIMU_BATCH simulates power system transients using the MATLAB Power System Toolbox
% This m-file takes the dynamic and load flow data and calculates the
% response of the power system to a fault (which is specified in a
% switching file) or other power system perturbance. See one of the supplied
% examples for file format and/or replacing technique.
%
%   NOTES:  Will run whatever is specified in the DataFile.m located in the
%           same folder as this script.
%
%   Runs scripts:
%   octaveComp      - Checks if running in Octave, takes compatibility steps
%   DataFile        - Contains simulation system information
%   handleNewGlobals - assigns system variables to global structure g
%   livePlot        - Plot data during simulation
%
%   Function calls to:
%   loadflow        - solve AC load flow
%   lfdcs           - solve load flow with DC lines
%   y_switch        - creates reduced Y matracies for fault scenarios
%
%   Create indicies in global g by calling:
%   mac_indx, exc_indx, tg_indx, dpwf_indx, pss_indx, svc_indx, tcsc_indx,
%   lm_indx, rlm_indx, pwrmod_indx, ...
%
%   Perform network and dynamic calculations by calling:
%   ...
%
%   Input:
%   N/A
%
%   Output:
%   Under revision...
%
%   History:
%   Date        Time    Engineer        Description
%   02/xx/97    XX:XX   Graham Rogers  	Version 1.0
%   08/xx/97    xx:xx   Graham Rogers   Version 1.1 - add modulation for load, exciter, and governor.
%   08/xx/97    XX:XX   Graham Rogers   Version 1.2 - add induction generator
%   (c) Copyright: Joe Chow/ Cherry Tree Scientific Software 1991 to 1997 - All rights reserved
%   06/xx/98    xx:xx   Graham Rogers   Version 1.3 - Add hydraulic turbine/governor
%   07/xx/98    XX:XX   Graham Rogers   Version 1.4 - add deltaP/omega filter
%   12/xx/98    xx:xx   Graham Rogers   Version 1.5 - add svc user defined damping controls.
%   01/xx/99    XX:XX   Graham Rogers   Version 1.5 - add tcsc model and tcsc user defined damping controls.
%   06/14/99    09:59   Graham Rogers   Version 1.6 - add user defined hvdc control
%                                       modify dc so that it is integrated with a sub-multiple time constant
%   09/xx/99    XX:XX   Graham Rogers   Version 1.7 - Add simple exciter with pi avr - smppi
%   xx/xx/15    xx:xx   Dan Trudnowski  Version 1.75 - Added pwrmod code
%   07/19/17    xx:xx   Rui             Version 1.8 - Add code for multiple HVDC systems
%   xx/xx/19    xx:xx   Dan Trudnowski  Version 1.9 - Added ivmmod code
%   07/15/20    14:27   Thad Haines     Version 2.0 - Revised format of globals and internal function documentation
%   07/16/20    11:19   Thad Haines     V 2.0.1 - Added cleanZeros to end of script to clean global g
%   07/18/20    10:52   Thad Haines     V 2.0.2 - Added Line monitor, area, and sytem average frequency calculations.
%   07/21/20    16:20   Thad haines     V 2.0.3 - Added AGC
%   07/29/20    15:20   Thad Haines     jay -> 1j

%%
%clear all
%clear global
% the above clears were removed to allow for running w/o running DataFile.m 5/20/20
% assumes required arrays created before this script runs and DataFile is delted/not found
format compact;
disp('*** s_simu_Batch Start')

close % close graphics windows
tic % set timer

disp('*** Declare Global Variables')
%% Contents of pst_var copied into this section so that globals highlight
% globals converted to the global g have been removed

% old method for declaring globals.
%pst_var % set up global variables (very many)

%% Remaining 'loose' globals

%% ivm variables - 5
global n_ivm mac_ivm_idx ivmmod_data ivmmod_d_sig ivmmod_e_sig

%% DeltaP/omega filter variables - 21
global  dpw_con dpw_out dpw_pot dpw_pss_idx dpw_mb_idx dpw_idx n_dpw dpw_Td_idx dpw_Tz_idx
global  sdpw1 sdpw2 sdpw3 sdpw4 sdpw5 sdpw6
global  dsdpw1 dsdpw2 dsdpw3 dsdpw4 dsdpw5 dsdpw6

%% pss design - 3 - Not used in Simulation? - thad 07/18/20
global ibus_con  netg_con  stab_con

%% global structured array
global g

%% unused/unimplemented damping controls -thad 07/15/20
% intentionally removed/ignored?
g.svc.svc_dc=[];

g.tcsc.tcsc_dc=[];
g.tcsc.n_tcscud = 0;

g.dc.dcr_dc=[];
g.dc.dci_dc=[];

%% input data file from d_*.m file
% 05/20 Edits - thad
% Check for Octave, automatically load compatibility script
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
g.mac.ibus_con = []; % ignore infinite buses in transient simulation 

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
    [bus_sol,line,~] = loadflow(g.bus.busOG,g.line.lineOG,tol,iter_max,acc,'n',2);
    g.bus.bus = bus_sol;  % solved loadflow solution needed for initialization
    g.line.line = line;
    clear bus_sol line
    %save sim_fle.mat bus line % no need in batch runs - thad 07/17/20
else
    % Has HVDC, use DC load flow
    [bus_sol,line,~,rec_par, inv_par, line_par] = lfdcs(g.bus.busOG,g.line.lineOG,g.dc.dci_dc,g.dc.dcr_dc);
    g.bus.bus = bus_sol;
    g.line.line = line;
    g.dc.rec_par = rec_par;
    g.dc.inv_par = inv_par;
    g.dc.line_par = line_par;
    clear bus_sol line rec_par inv_par line_par
    %save sim_fle.mat bus line rec_par  inv_par line_par% no need in batch runs - thad 07/17/20
end

%% set indexes
% note: dc index set in dc load flow
mac_indx();
exc_indx();
tg_indx();
dpwf_indx();
pss_indx();
svc_indx();
tcsc_indx();
lm_indx;
rlm_indx();
pwrmod_indx(g.bus.bus);
lmon_indx;
area_indx;
agc_indx;

% Handled in mac_indx
% g.ind.n_mot = size(g.ind.ind_con,1); % inductive motors
% g.igen.n_ig = size(g.igen.igen_con,1); % inductive generators
%
% if isempty(g.ind.n_mot)
%     g.ind.n_mot = 0;
% end
% if isempty(g.igen.n_ig)
%     g.igen.n_ig = 0;
% end

% ntot = g.mac.n_mac+g.ind.n_mot+g.igen.n_ig; % unused? commented out - thad 07/16/20
% ngm = g.mac.n_mac + g.ind.n_mot;

g.mac.n_pm = g.mac.n_mac; % used for pm modulation -- put into mac or tg indx?

%% Make sure bus max/min Q is the same as the pwrmod_con max/min Q - moved to place after pwrmod is counted (never actually executed prior...) -thad 06/30/20
if ~isempty(g.pwr.n_pwrmod)
    for kk=1:g.pwr.n_pwrmod
        n = find(g.pwr.pwrmod_con(kk,1)==g.bus.bus(:,1));
        g.bus.bus(n,11:12) = g.pwr.pwrmod_con(kk,6:7);
    end
    clear kk n
end

%% construct simulation switching sequence as defined in sw_con
warning('*** Initialize time and switching variables')

%tswitch(1) = g.sys.sw_con(1,1); -unused -thad 07/16/20

k = 1;
kdc = 1;

n_switch = length(g.sys.sw_con(:,1));

g.k.k_inc = zeros(n_switch-1,1);
g.k.k_incdc = zeros(n_switch-1,1);

t_switch = zeros(n_switch,1);
g.k.h = zeros(n_switch,1);
g.k.h_dc = zeros(n_switch,1);

for sw_count = 1:n_switch-1
    g.k.h(sw_count) = g.sys.sw_con(sw_count,7);%specified time step
    
    if g.k.h(sw_count)==0
        g.k.h(sw_count) = 0.01;
    end % default time step
    
    g.k.k_inc(sw_count) = fix((g.sys.sw_con(sw_count+1,1)-g.sys.sw_con(sw_count,1))/g.k.h(sw_count));%nearest lower integer
    
    if g.k.k_inc(sw_count)==0
        g.k.k_inc(sw_count)=1;
    end% minimum 1
    
    g.k.h(sw_count) = (g.sys.sw_con(sw_count+1,1)-g.sys.sw_con(sw_count,1))/g.k.k_inc(sw_count);%step length
    g.k.h_dc(sw_count) = g.k.h(sw_count)/10;
    g.k.k_incdc(sw_count) = 10*g.k.k_inc(sw_count);
    t_switch(sw_count+1) = t_switch(sw_count) +  g.k.k_inc(sw_count)*g.k.h(sw_count);
    t(k:k-1+g.k.k_inc(sw_count)) = t_switch(sw_count):g.k.h(sw_count):t_switch(sw_count+1)-g.k.h(sw_count);
    if ~isempty(g.dc.dcl_con)
        g.dc.t_dc(kdc:kdc-1+g.k.k_incdc(sw_count)) = t_switch(sw_count):g.k.h_dc(sw_count):t_switch(sw_count+1)-g.k.h_dc(sw_count);
    end
    k = k + g.k.k_inc(sw_count);
    kdc = kdc + g.k.k_incdc(sw_count);
end

% time for dc - multi-rate...
if ~isempty(g.dc.dcsp_con)
    g.dc.t_dc(kdc)=g.dc.t_dc(kdc-1)+g.k.h_dc(sw_count);
    for kk=1:10
        kdc=kdc+1;
        g.dc.t_dc(kdc)=g.dc.t_dc(kdc-1)+g.k.h_dc(sw_count);
    end
end
k = sum(g.k.k_inc)+1; % k is the total number of time steps in the simulation

t(k) = g.sys.sw_con(n_switch,1);

% add time array t to global g - thad
g.sys.t = t;

%% =====================================================================================================   
%% Start of Initializing Zeros =========================================================================
 
%% create zero matrices for variables to make algorithm more efficient?
warning('*** Initialize zero matricies...')

n = size(g.mac.mac_con, 1) ;
z = zeros(n,k);

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
        g.dc.Vdc(ire,:) = g.dc.rec_par(ihvdc_count,2);
        g.dc.i_dc(ire,:) = g.dc.line_par(ihvdc_count);    % for PDCI
        g.dc.i_dcr(ihvdc_count,:) = g.dc.i_dc(ire,:);
        g.dc.alpha(ihvdc_count,:) = g.dc.rec_par(ihvdc_count,1)*pi/180;
    end
    
    for ihvdc_count=1:g.dc.n_dcl
        iin=g.dc.i_idx(ihvdc_count);
        g.dc.Vdc(iin,:)= g.dc.inv_par(ihvdc_count,2);
        g.dc.i_dc(iin,:) = g.dc.line_par(ihvdc_count);
        g.dc.i_dci(ihvdc_count,:) = g.dc.i_dc(iin,:);
        g.dc.gamma(ihvdc_count,:) = g.dc.inv_par(ihvdc_count,1)*pi/180;
    end
    % end modification by Rui
    clear ihvdc_count
    
    g.dc.Vdc(g.dc.r_idx,:) = g.dc.rec_par(:,2);
    g.dc.Vdc(g.dc.i_idx,:) = g.dc.inv_par(:,2);
    g.dc.i_dc(g.dc.r_idx,:) = g.dc.line_par;
    g.dc.i_dc(g.dc.i_idx,:) = g.dc.line_par;
    g.dc.i_dcr(:,:) = g.dc.i_dc(g.dc.r_idx,:);
    g.dc.i_dci(:,:) = g.dc.i_dc(g.dc.i_idx,:);
    g.dc.alpha(:,:) = g.dc.rec_par(:,1)*pi/180;
    g.dc.gamma(:,:) = g.dc.inv_par(:,1)*pi/180;
    
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
clear j sv

%% End of DC specific stuff? - continuation of zero intialization - 06/01/20 -thad

%mac_ref = z1;  % unsure of this use
%sys_ref = z1;   % unsure of this use - thad 07/02/20

n_bus = length(g.bus.bus(:,1));
g.bus.bus_v = zeros(n_bus+1,k);

g.mac.cur_re = z;
g.mac.cur_im = z;
g.mac.psi_re = z;
g.mac.psi_im = z;

g.bus.theta = zeros(n_bus+1,k);

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

clear z_pss

%% delta P omwga filter states and derivative place holders -thad 07/15/20
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
clear z_dpw

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
        clear j
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
    g.svc.xsvc_dc = zeros(1,k);
    g.svc.dxsvc_dc = zeros(1,k);
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
        clear j
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

%% Initialize zeros for line monitor
if g.lmon.n_lmon ~=0
    for Lndx = 1:g.lmon.n_lmon
        g.lmon.line(Lndx).iFrom = zeros(1,k);
        g.lmon.line(Lndx).iTo = zeros(1,k);
        g.lmon.line(Lndx).sFrom = zeros(1,k);
        g.lmon.line(Lndx).sTo = zeros(1,k);
    end
end
clear Lndx

%% Initialize zeros for area values
if g.area.n_area ~=0
    for areaN = 1:g.area.n_area
        g.area.area(areaN).totH = zeros(1,k); % to account for possible trips
        g.area.area(areaN).aveF = zeros(1,k);
        g.area.area(areaN).totGen = zeros(1,k);
        g.area.area(areaN).totLoad = zeros(1,k);
        g.area.area(areaN).icA = zeros(1,k); % Actual interchange
        g.area.area(areaN).icS = ones(1,k); % Scheduled interchange
    end
end
clear areaN

%% Initialize zeros agc
if g.agc.n_agc ~=0
    for ndx = 1:g.agc.n_agc
        g.agc.agc(ndx).race = zeros(1,k);
        g.agc.agc(ndx).d_sace = zeros(1,k); % derivative
        g.agc.agc(ndx).sace = zeros(1,k);
        g.agc.agc(ndx).ace2dist = zeros(1,k);
        %g.agc.agc(ndx).iace = zeros(1,k); % window integrator not implemented
        
        for gndx=1:g.agc.agc(ndx).n_ctrlGen
            g.agc.agc(ndx).ctrlGen(gndx).input = zeros(1,k);
            g.agc.agc(ndx).ctrlGen(gndx).dx = zeros(1,k);
            g.agc.agc(ndx).ctrlGen(gndx).x = zeros(1,k);
        end
    end
end
clear ndx

%% Initialize zeros for logged system values
g.sys.aveF = zeros(1,k);
g.sys.totH = zeros(1,k);

 %% =====================================================================================================   
 %% Start Initialization ================================================================================
 

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


%% Flag = 0 == Initialization
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
if n_ivm~=0
    disp('ivm modulation')
    [~,~,~,~,Dini,Eini] = ivmmod_dyn([],[],g.bus.bus,g.sys.t,1,flag);
    
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
    
    g.k.h_sol = i_simu(1,1,g.k.k_inc,g.k.h,g.bus.bus_sim,Y1,Y2,Y3,Y4,Vr1,Vr2,bo);
    % reinitialize dc controls
    mdc_sig(1);
    dc_cont(0,1,1,g.bus.bus,flag);
    % initialize dc line
    dc_line(0,1,1,g.bus.bus,flag);
end

%H_sum = sum(g.mac.mac_con(:,16)./g.mac.mac_pot(:,1)); % unused?

%% step 3: Beginning of Huen's  (predictor-corrector) method
% Create indicies for simulation
kt = 0;
g.k.ks = 1;

k_tot = sum(g.k.k_inc);
lswitch = length(g.k.k_inc);
ktmax = k_tot-g.k.k_inc(lswitch);
g.bus.bus_sim = g.bus.bus;

% added from v2.3 06/01/20 - thad
g.mac.mac_trip_flags = false(g.mac.n_mac,1);
g.mac.mac_trip_states = 0;

%% Simulation loop start
warning('*** Simulation Loop Start')
while (kt<=ktmax)
    k_start = kt+1;
    
    if kt==ktmax
        k_end = kt + g.k.k_inc(g.k.ks);
    else
        k_end = kt + g.k.k_inc(g.k.ks) + 1;
    end
    
    for k = k_start:k_end

        % display k and t at g.k.k_inc and every ...th step - thad
        if ( mod(k,50)==0 ) || k == 1 || k == k_end
            fprintf('*** k = %5d, \tt(k) = %7.4f\n',k,g.sys.t(k)) % DEBUG
        end
        
     
        %% step 3a: network solution
         
        % mach_ref(k) = mac_ang(syn_ref,k);
        g.sys.mach_ref(k) = 0;
        g.mac.pmech(:,k+1) = g.mac.pmech(:,k);
        g.igen.tmig(:,k+1) = g.igen.tmig(:,k);
        
        if g.dc.n_conv~=0
            g.dc.cur_ord(:,k+1) = g.dc.cur_ord(:,k);
        end
        
        % Trip gen - Copied from v2.3 06/01/20 - thad
        [f,g.mac.mac_trip_states] = mac_trip_logic(g.mac.mac_trip_flags,g.mac.mac_trip_states,g.sys.t,k);
        g.mac.mac_trip_flags = g.mac.mac_trip_flags | f;
        
 %% =====================================================================================================   
 %% Start of Network Solution 1 =========================================================================
  
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
        
        %% Calculate current injections and bus voltages and angles
        if k >= sum(g.k.k_inc(1:3))+1
            %% fault cleared - post fault 2
            g.line.line_sim = g.line.line_pf2;
            g.bus.bus_sim = g.bus.bus_pf2;
            g.bus.bus_int = g.bus.bus_intpf2;
            
            Y1 = g.int.Y_gpf2;
            Y2 = g.int.Y_gncpf2;
            Y3 = g.int.Y_ncgpf2;
            Y4 = g.int.Y_ncpf2;
            Vr1 = g.int.V_rgpf2;
            Vr2 = g.int.V_rncpf2;
            bo = g.int.bopf2;
            
            % i_simu forms the network interface variables
            %h_sol = i_simu(k,ks,k_inc,h,bus_sim,Y1,Y2,Y3,Y4,Vr1,Vr2,bo);
            % duplicate call?
            % h_sol calculated after this 'if' block...
            
        elseif k >=sum(g.k.k_inc(1:2))+1
            %% near bus cleared - post fault 1
            g.line.line_sim = g.line.line_pf1;
            
            g.bus.bus_sim = g.bus.bus_pf1;
            g.bus.bus_int = g.bus.bus_intpf1;
            
            Y1 = g.int.Y_gpf1;
            Y2 = g.int.Y_gncpf1;
            Y3 = g.int.Y_ncgpf1;
            Y4 = g.int.Y_ncpf1;
            Vr1 = g.int.V_rgpf1;
            Vr2 = g.int.V_rncpf1;
            bo = g.int.bopf1;
            
            %h_sol = i_simu(k,ks,k_inc,h,bus_sim,Y1,Y2,Y3,Y4,Vr1,Vr2,bo);
            
        elseif k>=g.k.k_inc(1)+1
            %% fault applied - fault
            g.line.line_sim = g.line.line_f;
            g.bus.bus_sim = g.bus.bus_f;
            
            g.bus.bus_int = g.bus.bus_intf;
            
            Y1 = g.int.Y_gf;
            Y2 = g.int.Y_gncf;
            Y3 = g.int.Y_ncgf;
            Y4 = g.int.Y_ncf;
            Vr1 = g.int.V_rgf;
            Vr2 = g.int.V_rncf;
            bo = g.int.bof;
            
            %h_sol = i_simu(k,ks,k_inc,h,bus_sim,Y1,Y2,Y3,Y4,Vr1,Vr2,bo);
            
        elseif k<g.k.k_inc(1)+1
            %% pre fault
            g.line.line_sim = g.line.line;
            g.bus.bus_sim = g.bus.bus;
            
            g.bus.bus_int = g.bus.bus_intprf;
            
            Y1 = g.int.Y_gprf;
            Y2 = g.int.Y_gncprf;
            Y3 = g.int.Y_ncgprf;
            Y4 = g.int.Y_ncprf;
            Vr1 = g.int.V_rgprf;
            Vr2 = g.int.V_rncprf;
            bo = g.int.boprf;
            
            %h_sol = i_simu(k,ks,k_inc,h,bus_sim,Y1,Y2,Y3,Y4,Vr1,Vr2,bo);
        end
        
        %% apply gen trip - added from v2.3 - 06/01/20 - thad
        if sum(g.mac.mac_trip_flags)>0.5
            genBuses = g.mac.mac_con(g.mac.mac_trip_flags==1,2);
            for kB=1:length(genBuses)
                nL = find(genBuses(kB)==g.line.line_sim(:,1) | genBuses(kB)==g.line.line_sim(:,2));
                if isempty(nL); error(' '); end
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
        
 %% =====================================================================================================   
 %% Start Dynamic Solution 1 ============================================================================
 
        %% step 3b: compute dynamics and integrate
        flag = 2;
        %g.sys.sys_freq(k) = 1.0; % why?... 5/21/20 -thad
        % initialized as all ones on line 764
        
        mpm_sig(k);
        
        mac_ind(0,k,g.bus.bus_sim,flag);
        mac_igen(0,k,g.bus.bus_sim,flag);
        
        mac_sub(0,k,g.bus.bus_sim,flag);
        mac_tra(0,k,g.bus.bus_sim,flag);
        mac_em(0,k,g.bus.bus_sim,flag);
        
        dpwf(0,k,flag);
        pss(0,k,flag);
        
        mexc_sig(k);  % executed twice? - thad 07/18/20
        smpexc(0,k,flag);
        smppi(0,k,flag);
        exc_st3(0,k,flag);
        exc_dc12(0,k,flag);
        
        mtg_sig(k); % executed twice? - thad 07/18/20
        % AGC calculations after tg_sig, before tg dynamics
        if g.agc.n_agc ~= 0
            agc(k,flag); % perform calculations
        end
        
        tg(0,k,flag);
        tg_hydro(0,k,g.bus.bus_sim,flag);
        
        if g.svc.n_svc~=0
            v_svc = abs(g.bus.bus_v(g.bus.bus_int(g.svc.svc_con(:,2)),k));
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
            svc(0,k,g.bus.bus_sim,flag,v_svc);
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
            [~,~,dp,dq,~,~] = pwrmod_dyn(Pst,Qst,g.bus.bus,g.sys.t,k,flag,g.pwr.n_pwrmod);
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
            [P,Q,~,~] = pwrmod_dyn(Pst,Qst,g.bus.bus,g.sys.t,k,1,g.pwr.n_pwrmod); %update pwrmod_p_sig and pwrmod_q_sig
            if (length(P)~=g.pwr.n_pwrmod) || (length(Q)~=g.pwr.n_pwrmod)
                error('Dimension error in pwrmod_dyn');
            end
            g.pwr.pwrmod_p_sig(:,k) = P;
            g.pwr.pwrmod_q_sig(:,k) = Q;
            pwrmod_p(0,k,g.bus.bus_sim,flag);
            pwrmod_q(0,k,g.bus.bus_sim,flag);
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
            [d,e,~,~,~,~] = ivmmod_dyn(dst,est,g.bus.bus,g.sys.t,k,1); %get internal voltage signals
            if (length(d)~=n_ivm) || (length(e)~=n_ivm)
                error('Dimension error in ivmmod_dyn');
            end
            ivmmod_d_sig(:,k) = d;
            ivmmod_e_sig(:,k) = e;
            mac_ivm(0,k,g.bus.bus_sim,flag);
            [~,~,dd,de,~,~] = ivmmod_dyn(dst,est,g.bus.bus,g.sys.t,k,flag);
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
        
 %% =====================================================================================================   
 %% Start of DC solution 1 ==============================================================================
 
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
 
 %% =====================================================================================================   
 %% Start of Pretictor Integration ======================================================================
 
        %% following statements are predictor steps
        j = k+1;
        g.mac.mac_ang(:,j) = g.mac.mac_ang(:,k) + g.k.h_sol*g.mac.dmac_ang(:,k);
        g.mac.mac_spd(:,j) = g.mac.mac_spd(:,k) + g.k.h_sol*g.mac.dmac_spd(:,k);
        g.mac.edprime(:,j) = g.mac.edprime(:,k) + g.k.h_sol*g.mac.dedprime(:,k);
        g.mac.eqprime(:,j) = g.mac.eqprime(:,k) + g.k.h_sol*g.mac.deqprime(:,k);
        g.mac.psikd(:,j) = g.mac.psikd(:,k) + g.k.h_sol*g.mac.dpsikd(:,k);
        g.mac.psikq(:,j) = g.mac.psikq(:,k) + g.k.h_sol*g.mac.dpsikq(:,k);
        
        % Exciter integration
        g.exc.Efd(:,j) = g.exc.Efd(:,k) + g.k.h_sol*g.exc.dEfd(:,k);
        g.exc.V_R(:,j) = g.exc.V_R(:,k) + g.k.h_sol*g.exc.dV_R(:,k);
        g.exc.V_As(:,j) = g.exc.V_As(:,k) + g.k.h_sol*g.exc.dV_As(:,k);
        g.exc.R_f(:,j) = g.exc.R_f(:,k) + g.k.h_sol*g.exc.dR_f(:,k);
        g.exc.V_TR(:,j) = g.exc.V_TR(:,k) + g.k.h_sol*g.exc.dV_TR(:,k);
        
        if n_dpw ~= 0
            % only calculate if dpw filter is used
            sdpw1(:,j) = sdpw1(:,k) + g.k.h_sol*dsdpw1(:,k);
            sdpw2(:,j) = sdpw2(:,k) + g.k.h_sol*dsdpw2(:,k);
            sdpw3(:,j) = sdpw3(:,k) + g.k.h_sol*dsdpw3(:,k);
            sdpw4(:,j) = sdpw4(:,k) + g.k.h_sol*dsdpw4(:,k);
            sdpw5(:,j) = sdpw5(:,k) + g.k.h_sol*dsdpw5(:,k);
            sdpw6(:,j) = sdpw6(:,k) + g.k.h_sol*dsdpw6(:,k);
        end
        
        g.pss.pss1(:,j) = g.pss.pss1(:,k) + g.k.h_sol*g.pss.dpss1(:,k);
        g.pss.pss2(:,j) = g.pss.pss2(:,k) + g.k.h_sol*g.pss.dpss2(:,k);
        g.pss.pss3(:,j) = g.pss.pss3(:,k) + g.k.h_sol*g.pss.dpss3(:,k);
        
        % modified to g - thad
        g.tg.tg1(:,j) = g.tg.tg1(:,k) + g.k.h_sol*g.tg.dtg1(:,k);
        g.tg.tg2(:,j) = g.tg.tg2(:,k) + g.k.h_sol*g.tg.dtg2(:,k);
        g.tg.tg3(:,j) = g.tg.tg3(:,k) + g.k.h_sol*g.tg.dtg3(:,k);
        g.tg.tg4(:,j) = g.tg.tg4(:,k) + g.k.h_sol*g.tg.dtg4(:,k);
        g.tg.tg5(:,j) = g.tg.tg5(:,k) + g.k.h_sol*g.tg.dtg5(:,k);
        
        % induction motor integrations
        if g.ind.n_mot ~= 0
            g.ind.vdp(:,j) = g.ind.vdp(:,k) + g.k.h_sol*g.ind.dvdp(:,k);
            g.ind.vqp(:,j) = g.ind.vqp(:,k) + g.k.h_sol*g.ind.dvqp(:,k);
            g.ind.slip(:,j) = g.ind.slip(:,k) + g.k.h_sol*g.ind.dslip(:,k);
        end
        
        % induction generator integrations
        if g.igen.n_ig ~=0
            g.igen.vdpig(:,j) = g.igen.vdpig(:,k) + g.k.h_sol*g.igen.dvdpig(:,k);
            g.igen.vqpig(:,j) = g.igen.vqpig(:,k) + g.k.h_sol*g.igen.dvqpig(:,k);
            g.igen.slig(:,j) = g.igen.slig(:,k) + g.k.h_sol*g.igen.dslig(:,k);
        end
        
        % svc
        if g.svc.n_svc ~= 0
            g.svc.B_cv(:,j) = g.svc.B_cv(:,k) + g.k.h_sol*g.svc.dB_cv(:,k);
            g.svc.B_con(:,j) = g.svc.B_con(:,k) + g.k.h_sol*g.svc.dB_con(:,k);
            g.svc.xsvc_dc(:,j) = g.svc.xsvc_dc(:,k) + g.k.h_sol* g.svc.dxsvc_dc(:,k);
        end
        
        %tcsc
        if g.tcsc.n_tcsc ~= 0
            g.tcsc.B_tcsc(:,j) = g.tcsc.B_tcsc(:,k) + g.k.h_sol*g.tcsc.dB_tcsc(:,k);
            g.tcsc.xtcsc_dc(:,j) = g.tcsc.xtcsc_dc(:,k) + g.k.h_sol* g.tcsc.dxtcsc_dc(:,k);
        end
        
        if g.lmod.n_lmod~=0
            g.lmod.lmod_st(:,j) = g.lmod.lmod_st(:,k) + g.k.h_sol*g.lmod.dlmod_st(:,k); % line using g
        end
        
        if g.rlmod.n_rlmod~=0
            g.rlmod.rlmod_st(:,j) = g.rlmod.rlmod_st(:,k)+g.k.h_sol*g.rlmod.drlmod_st(:,k);
        end
        %% Copied from v2.3 - 06/01/20 - thad
        g.pwr.pwrmod_p_st(:,j) = g.pwr.pwrmod_p_st(:,k)+g.k.h_sol*g.pwr.dpwrmod_p_st(:,k);
        g.pwr.pwrmod_q_st(:,j) = g.pwr.pwrmod_q_st(:,k)+g.k.h_sol*g.pwr.dpwrmod_q_st(:,k);
        %% pwrmod
        if g.pwr.n_pwrmod~=0
            for index=1:g.pwr.n_pwrmod
                pwrmod_p_sigst{index}(:,j) = pwrmod_p_sigst{index}(:,k)+g.k.h_sol*dpwrmod_p_sigst{index}(:,k);
                pwrmod_q_sigst{index}(:,j) = pwrmod_q_sigst{index}(:,k)+g.k.h_sol*dpwrmod_q_sigst{index}(:,k);
            end
        end
        %% ivmmod
        if n_ivm~=0
            for index=1:n_ivm
                ivmmod_d_sigst{index}(:,j) = ivmmod_d_sigst{index}(:,k)+g.k.h_sol*divmmod_d_sigst{index}(:,k);
                ivmmod_e_sigst{index}(:,j) = ivmmod_e_sigst{index}(:,k)+g.k.h_sol*divmmod_e_sigst{index}(:,k);
            end
        end
        
        %% agc predictor integration
        if g.agc.n_agc ~=0
            for ndx = 1:g.agc.n_agc
                g.agc.agc(ndx).sace(j) = g.agc.agc(ndx).sace(k) + g.k.h_sol*g.agc.agc(ndx).d_sace(k);
                
                % integrate lowpass filter outs...
                for gndx=1:g.agc.agc(ndx).n_ctrlGen
                    g.agc.agc(ndx).ctrlGen(gndx).x(j) = g.agc.agc(ndx).ctrlGen(gndx).x(k)...
                        + g.k.h_sol * g.agc.agc(ndx).ctrlGen(gndx).dx(k)  ;
                end
            end
        end
        
 %% =====================================================================================================   
 %% Predictor Line Monitoring and Area Calculations =====================================================
         
        %% Line Monitoring
        if g.lmon.n_lmon~=0
            lmon(k)
        end
        
        %% Average Frequency Calculation
        calcAveF(k,1);
        
        %% Area Total Calcvulations
        if g.area.n_area ~= 0
            calcAreaVals(k,1);
        end
        
 %% =====================================================================================================   
 %% Start of Network Solution 2 =========================================================================
 
        %% Flag = 1
        % begining of solutions as j (i.e. Corrector step)
        flag = 1;
        % mach_ref(j) = mac_ang(syn_ref,j);
        g.sys.mach_ref(j) = 0;
        % perform network interface calculations again with predicted states
        mpm_sig(j);
        mac_ind(0,j,g.bus.bus_sim,flag);
        mac_igen(0,j,g.bus.bus_sim,flag);
        mac_sub(0,j,g.bus.bus_sim,flag);
        mac_tra(0,j,g.bus.bus_sim,flag);
        mac_em(0,j,g.bus.bus_sim,flag);
        mac_ivm(0,j,g.bus.bus_sim,flag);
        
        % assume Vdc remains unchanged for first pass through dc controls interface
        mdc_sig(j);
        dc_cont(0,j,10*(j-1)+1,g.bus.bus_sim,flag);
        
        % Calculate current injections and bus voltages and angles
        if j >= sum(g.k.k_inc(1:3))+1
            % fault cleared
            g.bus.bus_sim = g.bus.bus_pf2;
            g.bus.bus_int = g.bus.bus_intpf2;
            
            Y1 = g.int.Y_gpf2;
            Y2 = g.int.Y_gncpf2;
            Y3 = g.int.Y_ncgpf2;
            Y4 = g.int.Y_ncpf2;
            Vr1 = g.int.V_rgpf2;
            Vr2 = g.int.V_rncpf2;
            bo = g.int.bopf2;
            %h_sol = i_simu(j,ks,k_inc,h,bus_sim,Y1,Y2,Y3,Y4,Vr1,Vr2,bo);
            
        elseif j >=sum(g.k.k_inc(1:2))+1
            % near bus cleared
            g.bus.bus_sim = g.bus.bus_pf1;
            g.bus.bus_int = g.bus.bus_intpf1;
            
            Y1 = g.int.Y_gpf1;
            Y2 = g.int.Y_gncpf1;
            Y3 = g.int.Y_ncgpf1;
            Y4 = g.int.Y_ncpf1;
            Vr1 = g.int.V_rgpf1;
            Vr2 = g.int.V_rncpf1;
            bo = g.int.bopf1;
            %h_sol = i_simu(j,ks,k_inc,h,bus_sim,Y1,Y2,Y3,Y4,Vr1,Vr2,bo);
            
        elseif j>=g.k.k_inc(1)+1
            % fault applied
            g.bus.bus_sim = g.bus.bus_f;
            g.bus.bus_int = g.bus.bus_intf;
            
            Y1 = g.int.Y_gf;
            Y2 = g.int.Y_gncf;
            Y3 = g.int.Y_ncgf;
            Y4 = g.int.Y_ncf;
            Vr1 = g.int.V_rgf;
            Vr2 = g.int.V_rncf;
            bo = g.int.bof;
            %h_sol = i_simu(j,ks,k_inc,h,bus_sim,Y1,Y2,Y3,Y4,Vr1,Vr2,bo);
            
        elseif j<g.k.k_inc(1)+1  % JHC - DKF thinks k should be j.
            % Changed to j - thad 07/17/20
            % Though, since this is prefault - it should be handled by
            % previous elseif statements to handle the single corrector
            % step between a prefault and fault situation.
            
            % pre fault
            g.bus.bus_sim = g.bus.bus;
            g.bus.bus_int = g.bus.bus_intprf;
            
            Y1 = g.int.Y_gprf;
            Y2 = g.int.Y_gncprf;
            Y3 = g.int.Y_ncgprf;
            Y4 = g.int.Y_ncprf;
            Vr1 = g.int.V_rgprf;
            Vr2 = g.int.V_rncprf;
            bo = g.int.boprf;
            %h_sol = i_simu(j,ks,k_inc,h,bus_sim,Y1,Y2,Y3,Y4,Vr1,Vr2,bo);
        end
        
        % apply gen trip - copied from v2.3 - 06/01/20 - thad
        if sum(g.mac.mac_trip_flags)>0.5
            genBuses = g.mac.mac_con(g.mac.mac_trip_flags==1,2);
            for kB=1:length(genBuses)
                nL = find(genBuses(kB)==g.line.line_sim(:,1) | genBuses(kB)==g.line.line_sim(:,2));
                if isempty(nL)
                    error('nL is empty.');
                end
                g.line.line_sim(nL,4) = 1e7; %make reactance infinity
            end
            [Y1,Y2,Y3,Y4,Vr1,Vr2,bo] = red_ybus(g.bus.bus_sim,g.line.line_sim);
            clear nL kB genBuses
        end
        
        %% solve
        g.k.h_sol = i_simu(j,g.k.ks,g.k.k_inc,g.k.h,g.bus.bus_sim,Y1,Y2,Y3,Y4,Vr1,Vr2,bo);
        
        % only in step 2...
        g.mac.vex(:,j) = g.mac.vex(:,k);
        g.dc.cur_ord(:,j) = g.dc.cur_ord(:,k);
        
        % calculate the new value of bus angles rectifier user defined control
        if g.dc.ndcr_ud~=0
            tot_states=0;
            for jj = 1:g.dc.ndcr_ud
                b_num1 = g.dc.dcr_dc{jj,3};
                b_num2 = g.dc.dcr_dc{jj,4};
                conv_num = g.dc.dcr_dc{jj,2};
                g.dc.angdcr(jj,j) = g.bus.theta(g.bus.bus_int(b_num1),j)-g.bus.theta(g.bus.bus_int(b_num2),j);
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
                g.dc.angdci(jj,j) = g.bus.theta(g.bus.bus_int(b_num1),j)-g.bus.theta(g.bus.bus_int(b_num2),j);
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
        dc_cont(0,j,10*(j-1)+1,g.bus.bus_sim,flag);
        
        dpwf(0,j,flag);
        pss(0,j,flag);
        
        mexc_sig(j); % modulation
        smpexc(0,j,flag);
        smppi(0,j,flag);
        exc_st3(0,j,flag);
        exc_dc12(0,j,flag);
        
        tg(0,j,flag);
        tg_hydro(0,j,g.bus.bus_sim,flag);
        
        if g.svc.n_dcud~=0
            % set the new line currents
            for jj=1:g.svc.n_dcud
                l_num = g.svc.svc_dc{jj,3};
                svc_num = g.svc.svc_dc{jj,2};
                from_bus = g.bus.bus_int(g.line.line_sim(l_num,1));
                to_bus = g.bus.bus_int(g.line.line_sim(l_num,2));
                svc_bn = g.bus.bus_int(g.svc.svc_con(svc_num,2));
                V1 = g.bus.bus_v(from_bus,j);
                V2 = g.bus.bus_v(to_bus,j);
                R = g.line.line_sim(l_num,3);
                X = g.line.line_sim(l_num,4);
                B = g.line.line_sim(l_num,5);
                g.dc.tap = g.line.line_sim(l_num,6);
                phi = g.line.line_sim(l_num,7);
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
                g.tcsc.td_sig(jj,j) = abs(g.bus.bus_v(g.bus.bus_int(b_num),j));
            end
        end
        
 %% =====================================================================================================   
 %% Start of Dynamic Solution 2 =========================================================================
 
        %% Flag = 2, for 'corrector step' d's
        flag = 2;
        
        mpm_sig(j); % added for consistancy between steps
        
        mac_ind(0,j,g.bus.bus_sim,flag);
        mac_igen(0,j,g.bus.bus_sim,flag);
        
        mac_sub(0,j,g.bus.bus_sim,flag);
        mac_tra(0,j,g.bus.bus_sim,flag);
        mac_em(0,j,g.bus.bus_sim,flag);
        
        dpwf(0,j,flag);
        pss(0,j,flag);
        
        mexc_sig(j); % modulation
        smpexc(0,j,flag);
        smppi(0,j,flag);
        exc_st3(0,j,flag);
        exc_dc12(0,j,flag);
        
        mtg_sig(j);% modulation
        %% AGC calculations after tg_sig, before tg dynamics
        if g.agc.n_agc ~= 0
            agc(j,flag); % perform calculations
        end
        
        tg(0,j,flag);
        tg_hydro(0,j,g.bus.bus_sim,flag);
        
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
            v_svc = abs(g.bus.bus_v(g.bus.bus_int(g.svc.svc_con(:,2)),j));
            g.bus.bus_sim = svc(0,j,g.bus.bus_sim,flag,v_svc);
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
            ml_sig(j); % modulation
            lmod(0,j,flag);
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
            [~,~,dp,dq,~,~] = pwrmod_dyn(Pst,Qst,g.bus.bus,g.sys.t,j,flag,g.pwr.n_pwrmod);
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
            [P,Q,~,~,~,~] = pwrmod_dyn(Pst, Qst, g.bus.bus, g.sys.t, j, 1, g.pwr.n_pwrmod); %update pwrmod_p_sig and pwrmod_q_sig
            if (length(P)~=g.pwr.n_pwrmod) || (length(Q)~=g.pwr.n_pwrmod)
                error('Dimension error in pwrmod_dyn');
            end
            g.pwr.pwrmod_p_sig(:,j) = P;
            g.pwr.pwrmod_q_sig(:,j) = Q;
            pwrmod_p(0,j,g.bus.bus_sim,flag);
            pwrmod_q(0,j,g.bus.bus_sim,flag);
            clear P Q Pst Qst dp dq index
        end
        
        if n_ivm>0
            dst = cell(n_ivm,1);
            est = dst;
            for index=1:n_ivm
                dst{index} = ivmmod_d_sigst{index}(:,j);
                est{index} = ivmmod_e_sigst{index}(:,j);
            end
            [d,e,~,~,~,~] = ivmmod_dyn(dst,est,g.bus.bus,g.sys.t,j,1); % should this be g.bus.bus_sim? - thad 07/17/20
            if (length(d)~=n_ivm) || (length(e)~=n_ivm)
                error('Dimension error in ivmmod_dyn');
            end
            ivmmod_d_sig(:,j) = d;
            ivmmod_e_sig(:,j) = e;
            mac_ivm(0,j,g.bus.bus_sim,flag);
            [~,~,dd,de,~,~] = ivmmod_dyn(dst, est, g.bus.bus, g.sys.t, j, flag);
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
        
 %% =====================================================================================================   
 %% Start of DC Solution 2 ==============================================================================
                 
        %%integrate dc at ten times rate (DC Solution)
        if g.dc.n_conv~=0
            hdc_sol = g.k.h_sol/10;
            for kk = 1:10
                jdc=10*(j-1)+kk;
                [g.dc.xdcr_dc(:,jdc:jdc+1),g.dc.dxdcr_dc(:,jdc:jdc+1),g.dc.xdci_dc(:,jdc:jdc+1),g.dc.dxdci_dc(:,jdc:jdc+1)] = ...
                    dc_sim(j, kk, g.dc.dcr_dc, g.dc.dci_dc, g.dc.xdcr_dc(:,jdc), g.dc.xdci_dc(:,jdc), g.bus.bus_sim, hdc_sol);
            end
        else
            dc_cont(0, j, j, g.bus.bus_sim, 2);
            dc_line(0, j, j, g.bus.bus_sim, 2);
        end
        
 %% =====================================================================================================   
 %% Start of Corrector Integration ======================================================================
 
        %% following statements are corrector steps (RK2 computation)
        g.mac.mac_ang(:,j) = g.mac.mac_ang(:,k) + g.k.h_sol*(g.mac.dmac_ang(:,k)+g.mac.dmac_ang(:,j))/2.;
        g.mac.mac_spd(:,j) = g.mac.mac_spd(:,k) + g.k.h_sol*(g.mac.dmac_spd(:,k)+g.mac.dmac_spd(:,j))/2.;
        g.mac.edprime(:,j) = g.mac.edprime(:,k) + g.k.h_sol*(g.mac.dedprime(:,k)+g.mac.dedprime(:,j))/2.;
        g.mac.eqprime(:,j) = g.mac.eqprime(:,k) + g.k.h_sol*(g.mac.deqprime(:,k)+g.mac.deqprime(:,j))/2.;
        g.mac.psikd(:,j) = g.mac.psikd(:,k) + g.k.h_sol*(g.mac.dpsikd(:,k)+g.mac.dpsikd(:,j))/2.;
        g.mac.psikq(:,j) = g.mac.psikq(:,k) + g.k.h_sol*(g.mac.dpsikq(:,k)+g.mac.dpsikq(:,j))/2.;
        
        % exciter integration
        g.exc.Efd(:,j) = g.exc.Efd(:,k) + g.k.h_sol*(g.exc.dEfd(:,k)+g.exc.dEfd(:,j))/2.;
        g.exc.V_R(:,j) = g.exc.V_R(:,k) + g.k.h_sol*(g.exc.dV_R(:,k)+g.exc.dV_R(:,j))/2.;
        g.exc.V_As(:,j) = g.exc.V_As(:,k) + g.k.h_sol*(g.exc.dV_As(:,k)+g.exc.dV_As(:,j))/2.;
        g.exc.R_f(:,j) = g.exc.R_f(:,k) + g.k.h_sol*(g.exc.dR_f(:,k)+g.exc.dR_f(:,j))/2.;
        g.exc.V_TR(:,j) = g.exc.V_TR(:,k) + g.k.h_sol*(g.exc.dV_TR(:,k)+g.exc.dV_TR(:,j))/2.;
        
        % removed extra 1 in global names. - thad 07/06/20
        if n_dpw ~= 0
            % only calculate if dpw filter is used
            sdpw1(:,j) = sdpw1(:,k) +g.k.h_sol*(dsdpw1(:,k)+dsdpw1(:,j))/2.;
            sdpw2(:,j) = sdpw2(:,k) +g.k.h_sol*(dsdpw2(:,k)+dsdpw2(:,j))/2.;
            sdpw3(:,j) = sdpw3(:,k) +g.k.h_sol*(dsdpw3(:,k)+dsdpw3(:,j))/2.;
            sdpw4(:,j) = sdpw4(:,k) +g.k.h_sol*(dsdpw4(:,k)+dsdpw4(:,j))/2.;
            sdpw5(:,j) = sdpw5(:,k) +g.k.h_sol*(dsdpw5(:,k)+dsdpw5(:,j))/2.;
            sdpw6(:,j) = sdpw6(:,k) +g.k.h_sol*(dsdpw6(:,k)+dsdpw6(:,j))/2.;
        end
        
        g.pss.pss1(:,j) = g.pss.pss1(:,k) +g.k.h_sol*(g.pss.dpss1(:,k)+g.pss.dpss1(:,j))/2.;
        g.pss.pss2(:,j) = g.pss.pss2(:,k) +g.k.h_sol*(g.pss.dpss2(:,k)+g.pss.dpss2(:,j))/2.;
        g.pss.pss3(:,j) = g.pss.pss3(:,k) +g.k.h_sol*(g.pss.dpss3(:,k)+g.pss.dpss3(:,j))/2.;
        
        % modified to g
        g.tg.tg1(:,j) = g.tg.tg1(:,k) + g.k.h_sol*(g.tg.dtg1(:,k) + g.tg.dtg1(:,j))/2.;
        g.tg.tg2(:,j) = g.tg.tg2(:,k) + g.k.h_sol*(g.tg.dtg2(:,k) + g.tg.dtg2(:,j))/2.;
        g.tg.tg3(:,j) = g.tg.tg3(:,k) + g.k.h_sol*(g.tg.dtg3(:,k) + g.tg.dtg3(:,j))/2.;
        g.tg.tg4(:,j) = g.tg.tg4(:,k) + g.k.h_sol*(g.tg.dtg4(:,k) + g.tg.dtg4(:,j))/2.;
        g.tg.tg5(:,j) = g.tg.tg5(:,k) + g.k.h_sol*(g.tg.dtg5(:,k) + g.tg.dtg5(:,j))/2.;
        
        % induction motor integrations
        if g.ind.n_mot ~= 0
            g.ind.vdp(:,j) = g.ind.vdp(:,k) + g.k.h_sol*(g.ind.dvdp(:,j) + g.ind.dvdp(:,k))/2.;
            g.ind.vqp(:,j) = g.ind.vqp(:,k) + g.k.h_sol*(g.ind.dvqp(:,j) + g.ind.dvqp(:,k))/2.;
            g.ind.slip(:,j) = g.ind.slip(:,k) + g.k.h_sol*(g.ind.dslip(:,j) + g.ind.dslip(:,k))/2.;
        end
        
        % induction generator integrations
        if g.igen.n_ig ~=0
            g.igen.vdpig(:,j) = g.igen.vdpig(:,k) + g.k.h_sol*(g.igen.dvdpig(:,j) + g.igen.dvdpig(:,k))/2.;
            g.igen.vqpig(:,j) = g.igen.vqpig(:,k) + g.k.h_sol*(g.igen.dvqpig(:,j) + g.igen.dvqpig(:,k))/2.;
            g.igen.slig(:,j) = g.igen.slig(:,k) + g.k.h_sol*(g.igen.dslig(:,j) + g.igen.dslig(:,k))/2.;
        end
        
        % svc
        if g.svc.n_svc ~= 0
            g.svc.B_cv(:,j) = g.svc.B_cv(:,k) + g.k.h_sol*(g.svc.dB_cv(:,j) + g.svc.dB_cv(:,k))/2.;
            g.svc.B_con(:,j) = g.svc.B_con(:,k) + g.k.h_sol*(g.svc.dB_con(:,j) + g.svc.dB_con(:,k))/2.;
            g.svc.xsvc_dc(:,j) = g.svc.xsvc_dc(:,k) + g.k.h_sol*(g.svc.dxsvc_dc(:,j) + g.svc.dxsvc_dc(:,k))/2.;
        end
        
        %tcsc
        if g.tcsc.n_tcsc ~= 0
            g.tcsc.B_tcsc(:,j) = g.tcsc.B_tcsc(:,k) + g.k.h_sol*(g.tcsc.dB_tcsc(:,j) + g.tcsc.dB_tcsc(:,k))/2.;
            g.tcsc.xtcsc_dc(:,j) = g.tcsc.xtcsc_dc(:,k) + g.k.h_sol*(g.tcsc.dxtcsc_dc(:,j) + g.tcsc.dxtcsc_dc(:,k))/2.;
        end
        
        if g.lmod.n_lmod~=0
            g.lmod.lmod_st(:,j) = g.lmod.lmod_st(:,k) + g.k.h_sol*(g.lmod.dlmod_st(:,j) + g.lmod.dlmod_st(:,k))/2.; % modified line with g
        end
        if g.rlmod.n_rlmod~=0
            g.rlmod.rlmod_st(:,j) = g.rlmod.rlmod_st(:,k) + g.k.h_sol*(g.rlmod.drlmod_st(:,j) + g.rlmod.drlmod_st(:,k))/2.;
        end
        
        % Copied from v2.3 - 06/01/20 - thad
        g.pwr.pwrmod_p_st(:,j) = g.pwr.pwrmod_p_st(:,k)+g.k.h_sol*(g.pwr.dpwrmod_p_st(:,j) + g.pwr.dpwrmod_p_st(:,k))/2;
        g.pwr.pwrmod_q_st(:,j) = g.pwr.pwrmod_q_st(:,k)+g.k.h_sol*(g.pwr.dpwrmod_q_st(:,j) + g.pwr.dpwrmod_q_st(:,k))/2;
        if g.pwr.n_pwrmod~=0
            for index=1:g.pwr.n_pwrmod
                % make global? -thad 07/06/20
                pwrmod_p_sigst{index}(:,j) = pwrmod_p_sigst{index}(:,k)+g.k.h_sol*(dpwrmod_p_sigst{index}(:,j) + dpwrmod_p_sigst{index}(:,k))/2;
                pwrmod_q_sigst{index}(:,j) = pwrmod_q_sigst{index}(:,k)+g.k.h_sol*(dpwrmod_q_sigst{index}(:,j) + dpwrmod_q_sigst{index}(:,k))/2;
            end
        end
        
        if n_ivm~=0
            for index=1:n_ivm
                % make global? -thad 07/06/20
                ivmmod_d_sigst{index}(:,j) = ivmmod_d_sigst{index}(:,k)+g.k.h_sol*(divmmod_d_sigst{index}(:,j) + divmmod_d_sigst{index}(:,k))/2;
                ivmmod_e_sigst{index}(:,j) = ivmmod_e_sigst{index}(:,k)+g.k.h_sol*(divmmod_e_sigst{index}(:,j) + divmmod_e_sigst{index}(:,k))/2;
            end
        end
        
        %% agc corrector integration
        if g.agc.n_agc ~=0
            for ndx = 1:g.agc.n_agc
                g.agc.agc(ndx).sace(j) = g.agc.agc(ndx).sace(k) + g.k.h_sol*(g.agc.agc(ndx).d_sace(j)+g.agc.agc(ndx).d_sace(k))/2;
                
                % integrate lowpass filter outs...
                for gndx=1:g.agc.agc(ndx).n_ctrlGen
                    g.agc.agc(ndx).ctrlGen(gndx).x(j) = g.agc.agc(ndx).ctrlGen(gndx).x(k)...
                        + g.k.h_sol *(g.agc.agc(ndx).ctrlGen(gndx).dx(j)+  g.agc.agc(ndx).ctrlGen(gndx).dx(k))/2  ;
                end
            end
        end
        
 %% =====================================================================================================   
 %% Corrector Line Monitoring and Area Calculations =====================================================  
        
        %% Line Monitoring
        if g.lmon.n_lmon~=0
            lmon(j)
        end
        %% Average Frequency Calculation
        calcAveF(j,1);
        
        %% Area Total Calculations
        if g.area.n_area ~= 0
            calcAreaVals(j,1);
        end
        
 %% =====================================================================================================   
 %% end  of j step calculations =========================================================================
         
        %% Live plot call
        if g.sys.livePlotFlag
            livePlot(k)
        end
        

    end % end of k= k_start:k_end;
    % counter increment
    kt = kt + g.k.k_inc(g.k.ks);
    g.k.ks = g.k.ks+1;
end% end simulation loop

%% Final 'live' plot call
if g.sys.livePlotFlag
    livePlot('end')
end

% Now handled during simulation via lmon - below code can be used to verify functionality
% %% calculation of line currents post sim
% V1 = g.bus.bus_v(g.bus.bus_int(g.line.line(:,1)),:);
% V2 = g.bus.bus_v(g.bus.bus_int(g.line.line(:,2)),:);
% R = g.line.line(:,3);
% X = g.line.line(:,4);
% B = g.line.line(:,5);
% g.dc.tap = g.line.line(:,6);
% phi = g.line.line(:,7);
%
% [ilf,ilt] = line_cur(V1,V2,R,X,B,g.dc.tap,phi);%line currents
% [sInjF,sInjT] = line_pq(V1,V2,R,X,B,g.dc.tap,phi);% 'line flows' - complex power injection at bus

%% full sim timing
disp('*** End simulation.')
et = toc;
ets = num2str(et);
g.sys.ElapsedNonLinearTime = ets;
disp(['*** Elapsed Simulation Time = ' ets 's'])

%% Clean up logged DC variables to length of DC simulated time.
if ~isempty(g.dc.dcsp_con)
    disp('*** Adjusting logged data lengths...')
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
end

%% Clean up various temporary and function input values
clear V1 V2 R X B jj phi % used in calls to line_cur, line_pq
clear Y1 Y2 Y3 Y4 Vr1 Vr2 bo % used in calls to i_simu, red_ybus
clear et ets % used to store/display elapsed time

%% Remove all zero only data
varNames = who()'; % all variable names in workspace
clearedVars = {}; % cell to hold names of deleted 'all zero' variables

for vName = varNames
    try
        zeroTest = eval(sprintf('all(%s(:)==0)', vName{1})); % check if all zeros
        if zeroTest
            eval(sprintf('clear %s',vName{1}) ); % clear variable
            clearedVars{end+1} = vName{1}; % add name to cell for reference
        end
    catch %ME
        % catches stuctures
        testStr = ['isstruct(',vName{1},')'];
        if eval(testStr)
            fprintf('*** Clearing zeros from structure %s...\n', vName{1});
            eval(sprintf('[%s, clearedVars] = cleanZeros(%s, clearedVars);',vName{1},vName{1} ))
        end
        clear testStr
    end
    
end
g.sys.clearedVars = clearedVars; % attach cleard vars to global g
clear varNames vName zeroTest clearedVars % variables associated with clearing zeros.
%%
disp('*** s_simu_Batch End')
disp(' ')
