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
%   07/29/20    15:20   Thad Haines     jay -> 1j
%   08/11/20    11:36   Thad Haines     Added ivm to global

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

g.mac.n_pm = g.mac.n_mac; % used for pm modulation -- put into mac or tg indx?

%% Make sure bus max/min Q is the same as the pwrmod_con max/min Q
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
initZeros(k, kdc)

%% =====================================================================================================
%% Start Initialization ================================================================================
initNLsim()
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
        j = k+1;
        
        % display k and t at g.k.k_inc and every ...th step - thad
        if ( mod(k,50)==0 ) || k == 1 || k == k_end
            fprintf('*** k = %5d, \tt(k) = %7.4f\n',k,g.sys.t(k)) % DEBUG
        end
        
        %% Time step start
        initStep(k)
        
        %% Predictor Solution =============================================
        networkSolution(k)
        dynamicSolution(k)
        dcSolution(k)
        predictorIntegration(k, j, g.k.h_sol)
        monitorSolution(k);
        
        %% Corrector Solution =============================================
        networkSolution(j)
        dynamicSolution(j)
        dcSolution(j)      
        correctorIntegration(k, j, g.k.h_sol)
        monitorSolution(j);        
        
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
disp('*** s_simu_BatchTestF End')
disp(' ')
