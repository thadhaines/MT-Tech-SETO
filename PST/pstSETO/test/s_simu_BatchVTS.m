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
%   07/23/20    11:24   Thad Haines     Begining work on functionalization of solution.
%   07/27/20    09:54   Thad Haines     Begining work on Variable Time step simulation

%%
%clear all
%clear global
% the above clears were removed to allow for running w/o running DataFile.m 5/20/20
% assumes required arrays created before this script runs and DataFile is delted/not found
format compact;
disp('*** s_simu_BatchVTS Start')

close % close graphics windows
tic % set timer

jay = sqrt(-1);

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
% allocate more time by decreasing sw_con ts...
g.sys.sw_con(:,7) = g.sys.sw_con(:,7)./20;

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
initZeros(k, kdc)

%% Initialize required VTS globals
g.vts.dataN = 1;
handleStDx(1, 0, 0) % init

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
g.bus.bus_v(1:n_bus,1) = g.bus.bus(:,2).*exp(jay*g.bus.theta(1:n_bus,1)); %complex bus voltage

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
    i_acr = (Pr-jay*Qr)./conj(VLT(g.dc.r_idx));
    i_aci = (Pi - jay*Qi)./conj(VLT(g.dc.i_idx));
    IHT(g.dc.r_idx,1)=i_acr;
    IHT(g.dc.i_idx,1)=i_aci;
    g.dc.VHT(g.dc.r_idx,1) = (VLT(g.dc.r_idx) + jay*g.dc.dcc_pot(:,2).*i_acr);
    g.dc.VHT(g.dc.i_idx,1) = (VLT(g.dc.i_idx) + jay*g.dc.dcc_pot(:,4).*i_aci);
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
    
    g.k.h_sol = i_simu(1,1,g.k.k_inc, g.k.h, g.bus.bus_sim,Y1,Y2,Y3,Y4,Vr1,Vr2,bo);
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

%% creation of time blocks
initTblocks()
%% defining ODE input function
inputFcn = str2func('vtsInputFcn');
outputFcn = str2func('vtsOutputFcn');

% Configure ODE settings
%options = odeset('RelTol',1e-3,'AbsTol',1e-6); % default settings
%options = odeset('RelTol',1e-3,'AbsTol',1e-3, 'OutputFcn',outputFcn); % default settings
options = odeset('RelTol',1e-3,'AbsTol',1e-6, ... %'RelTol',1e-5,'AbsTol',1e-8, ...
    'InitialStep', 1/60/4, ...
    'MaxStep',20, ...
    'OutputFcn',outputFcn); % set 'OutputFcn' to function handle

handleStDx(1, [], 3) % update g.vts.stVec to initial conditions


%% Simulation loop start
warning('*** Simulation Loop Start')
g.vts.dataN = 1;
g.vts.iter = 0; % for keeping track of iteration...
g.vts.tot_iter = 0;
for simTblock = 1:size(g.vts.t_block)
    g.vts.t_blockN = simTblock;
    % 15s - slower during transients - faster when no action.
    % 113 - works well during transierts, slower during no action 
    % ode23s - many iterations per step - not very viable
    % ode23tb - occasionally hundereds of iterations, sometimes not... decent
    % ode23 - similar to 23tb, timstep doesn't get very large
    % ode23t - works...
    ode15s(inputFcn, g.vts.t_block(g.vts.t_blockN,:), g.vts.stVec , options);
    %feval(odeName, inputFcn, g.vts.t_block(simTblock,:), g.vts.stVec , options);
       
end% end simulation loop

% handle extraneous data index addition
g.vts.dataN = g.vts.dataN-1;

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
disp(['*** Total solutions = ' int2str(g.vts.tot_iter)])
disp(['*** Total data points = ' int2str(g.vts.dataN)])

%% TODO: trim logged values to length of g.vts.dataN

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
disp('*** s_simu_BatchTestF End')
disp(' ')
