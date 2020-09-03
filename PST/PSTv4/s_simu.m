%S_SIMU runs a non-linear simulation of a given data file.
% S_SIMU simulates power system transients using Power System Toolbox
% This m-file takes the dynamic and load flow data and calculates the
% response of the power system to a fault (which is specified in a
% switching file) or other power system perturbance. See one of the supplied
% examples for file format and/or replacing technique.
%
%   NOTES:  To run in stand-alone mode, clear all variables before
%           running script from PST main directory.
%           Accomodates for batch mode runs automatically.
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
%   ...
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
%   g - global structured variable containing trimmed and cleaned data from simulation.
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
%   07/16/20    11:19   Thad Haines     V 2.1 - Added cleanZeros to end of script to clean global g
%   07/18/20    10:52   Thad Haines     V 2.2 - Added Line monitor, area, and sytem average frequency calculations.
%   07/21/20    16:20   Thad Haines     V 2.3 - Added AGC
%   07/23/20    11:24   Thad Haines     V 2.3.1 - Begining work on functionalization of solution.
%   07/27/20    09:54   Thad Haines     V 2.3.2 - Begining work on Variable Time step simulation
%   07/28/20    16:20   Thad Haines     V 2.4 - Initial working VTS simulation
%   07/29/20    15:20   Thad Haines     V 2.4.1 - jay -> 1j
%   07/30/20    08:04   Thad Haines     V 2.5 - Re-integrated interactive script running/plotting
%   07/30/20    20:34   Thad Haines     V 2.6 - Added selectable solution methods
%   08/11/20    11:36   Thad Haines     V 2.7 - added ivm to global
%   08/13/20    09:37   Thad Haines     V 2.8 - Incorporated new license
%   08/16/20    20:28   Thad Haines     V 2.8.1 - Increase logged data only if VTS
%   08/18/20    10:51   Thad Haines     V 2.8.2 - Moved network solution in Huen's method to work with AGC
%   08/21/20    12:53   Thad Haines     V 2.8.3 - Handled FTS->VTS unique time vector issue
%   08/25/20    11:02   Thad Haines     V 2.8.4 - Added try catch to handle non-convergence
%   08/28/20    14:05   Thad Haines     V 2.8.5 - removed mac_trip_flags init, changed number of steps to print data

% (c) Montana Technological University / Thad Haines 2020
% (c) Montana Technological University / Daniel Trudnowski 2019
% (c) Montana Tech / Daniel Trudnowski 2015
% (c) Copyright: Joe Chow/ Cherry Tree Scientific Software 1991 to 2020 - All rights reserved
%
% Permission is hereby granted, free of charge, to any person obtaining a copy
% of this software and associated documentation files (the "Software"), to deal
% in the Software without restriction, including without limitation the rights
% to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
% copies of the Software, and to permit persons to whom the Software is
% furnished to do so.
%
% The above copyright notice and this permission notice shall be included in all
% copies or substantial portions of the Software.
%
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
% IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
% FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
% AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
% LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
% OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
% SOFTWARE.
%

format compact;
disp('***    PST v4.0.0-a.9    ***')
disp('*** s_simu Start')
disp('*** Declaring Global Variables...')

%% Remaining 'loose' globals
%% DeltaP/omega filter variables - 21 - model use not documented
global  dpw_con dpw_out dpw_pot dpw_pss_idx dpw_mb_idx dpw_idx n_dpw dpw_Td_idx dpw_Tz_idx
global  sdpw1 sdpw2 sdpw3 sdpw4 sdpw5 sdpw6
global  dsdpw1 dsdpw2 dsdpw3 dsdpw4 dsdpw5 dsdpw6

%% pss design - 3 - Not used in Simulation? - thad 07/18/20
global ibus_con  netg_con  stab_con

%% global structured array
global g

%% Check for Octave, automatically load compatibility script
% Assumes license of Octave will be 'GNU ...'
dataCheck = license;
if all(dataCheck(1:3)=='GNU')
    fprintf('*** Octave detected, loading compatiblity commands and packages...\n')
    octaveComp
else
    fprintf('*** MATLAB detected.\n')
end
clear dataCheck

%% IF run as a standalone script: querry user
% Assumes all variables up to this point are global.

if all( size(whos('global')) == size(whos()) )
    scriptRunFlag = 1;
    fprintf('*** Collecting system information from user...\n')
    % input data file
    [fname, pathname]=uigetfile('*.m','Select Data File to load');
    if ~ischar(fname) || ~ischar(pathname)
        error('Data File Not selected - Stopping simulation')
    end
    % replace data file in current directory with selected file.
    copyfile([pathname, fname],'DataFile.m'); % copy system data file to batch run location
    % ask for Fbase
    prompt={'System Frequency Base [Hz]:','System S Base [MW]:','Live Plot?' };
    name='Input F and S base values and optional Live Plot Flag';
    numlines=1;
    defaultanswer={'60','100','1'};
    
    answer=inputdlg(prompt,name,numlines,defaultanswer);
    
    Fbase = str2num(answer{1});
    Sbase = str2num(answer{2});
    livePlotFlag = str2num(answer{3});
    clear fname pathname prompt name numlines defaultanswer answer
else
    scriptRunFlag = 'No thanks - I like using batch scripts...';
end

%% account for non-existant DataFile
try
    DataFile %Batch name for data file
catch ME
    fprintf('*** Caught ERROR: %s\n',ME.message)
    fprintf('*** Continuing with simulation...\n')
end

% remove copied Data file
if scriptRunFlag == 1
    delete('DataFile.m')
end

%% run script to handle legacy input to new global g approach
handleNewGlobals

% check for valid dynamic data file
if isempty(g.mac.mac_con)
    error('! mac_con is Empty - invalid/incomplete input data!')
end
if isempty(g.sys.sw_con)
    error('! sw_con is Empty - simulation has no switching data!')
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
g.sys.basrad = 2*pi*g.sys.sys_freq;

% Handle variable input base MVA
% assumes Sbase defined in DataFile or earlier, basmva is defined as global in pst_var.
if ~exist('Sbase','var')
    fprintf('*** Sbase Not defined - assuming 100 MVA base.\n')
    g.sys.basmva = 100;
elseif isnumeric(Sbase)
    fprintf('*** Sbase found - Power base is set to %3.3f MVA\n', Sbase)
    g.sys.basmva = Sbase;
end

clear Fbase Sbase

%% other init operations
tic                     % start timer
%startStr = datestr(now);% start time String
%startClock = clock;     % start date time vector
g.sys.syn_ref = 0 ;     % synchronous reference frame
g.mac.ibus_con = [];    % ignore infinite buses in transient simulation

%% unused/unimplemented damping controls -thad 07/15/20
% intentionally removed/ignored?
g.svc.svc_dc = [];
g.tcsc.tcsc_dc = [];
g.tcsc.n_tcscud = 0;
g.dc.dcr_dc = [];
g.dc.dci_dc = [];

%% solve for loadflow - loadflow parameter
warning('*** Solve initial loadflow')
if isempty(g.dc.dcsp_con)
    % AC power flow
    g.dc.n_conv = 0;
    g.dc.n_dcl = 0;
    g.dc.ndcr_ud = 0;
    g.dc.ndci_ud = 0;
    tol = 1e-9;   % tolerance for convergence
    iter_max = 30; % maximum number of iterations
    acc = 1.0;   % acceleration factor
    [bus_sol,line,~] = loadflow(g.bus.busOG,g.line.lineOG,tol,iter_max,acc,'n',2);
    g.bus.bus = bus_sol;  % solved loadflow solution needed for initialization
    g.line.line = line;
    clear bus_sol line tol iter_max acc
else
    % Has HVDC, use DC load flow
    [bus_sol,line,~,rec_par, inv_par, line_par] = lfdcs(g.bus.busOG,g.line.lineOG,g.dc.dci_dc,g.dc.dcr_dc);
    g.bus.bus = bus_sol;
    g.line.line = line;
    g.dc.rec_par = rec_par;
    g.dc.inv_par = inv_par;
    g.dc.line_par = line_par;
    clear bus_sol line rec_par inv_par line_par
end

%% Set indexes
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
        g.bus.bus(n, 11:12) = g.pwr.pwrmod_con(kk,6:7);
    end
    clear kk n
end

%% construct simulation switching sequence as defined in sw_con (original method)
warning('*** Initialize time and switching variables')

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
    
    if g.k.h(sw_count) == 0
        g.k.h(sw_count) = 0.01;
    end % default time step
    
    % number of steps in 'time block'
    g.k.k_inc(sw_count) = fix((g.sys.sw_con(sw_count+1,1)-g.sys.sw_con(sw_count,1))/g.k.h(sw_count));%nearest lower integer
    
    if g.k.k_inc(sw_count) == 0
        g.k.k_inc(sw_count) = 1;
    end% minimum 1
    
    % adjust time step so integer number of steps in block
    g.k.h(sw_count) = (g.sys.sw_con(sw_count+1,1)-g.sys.sw_con(sw_count,1))/g.k.k_inc(sw_count);%step length
    g.k.h_dc(sw_count) = g.k.h(sw_count)/10; % h_dc 10 times faster
    g.k.k_incdc(sw_count) = 10*g.k.k_inc(sw_count);
    
    t_switch(sw_count+1) = t_switch(sw_count) +  g.k.k_inc(sw_count)*g.k.h(sw_count);
    
    % create time vector block from start time to end time - 1 timestep
    t(k:k-1+g.k.k_inc(sw_count)) = t_switch(sw_count):g.k.h(sw_count):t_switch(sw_count+1)-g.k.h(sw_count);
    
    if ~isempty(g.dc.dcl_con)
        g.dc.t_dc(kdc:kdc-1+g.k.k_incdc(sw_count)) = t_switch(sw_count):g.k.h_dc(sw_count):t_switch(sw_count+1)-g.k.h_dc(sw_count);
    end
    
    % keep track of indicies
    k = k + g.k.k_inc(sw_count);
    kdc = kdc + g.k.k_incdc(sw_count);
end

% time for dc - multi-rate
if ~isempty(g.dc.dcsp_con)
    g.dc.t_dc(kdc)=g.dc.t_dc(kdc-1)+g.k.h_dc(sw_count);
    for kk = 1:10
        kdc = kdc+1;
        g.dc.t_dc(kdc) = g.dc.t_dc(kdc-1)+g.k.h_dc(sw_count);
    end
    g.dc.t_dc_OLD = g.dc.t_dc;
end

k = sum(g.k.k_inc)+1; % k is the total number of time steps in the simulation (+1 for final predictor step)
t(k) = g.sys.sw_con(n_switch,1); % final time into time vector

% NOTE: the time vector created above is used only to create zeros...
g.sys.t = zeros(1,size(t,2));
g.sys.t_OLD = t; % for reference only

%% creation of VTS time blocks
initTblocks()

% If variable step, increase amount of zeros to log
if ~all(strcmp(g.vts.solver_con, 'huens'))
    k = max(size(g.sys.t_OLD))*4;
    g.sys.t = zeros(1,k);
    kdc = k*10;
else
    k = max(size(g.sys.t))+1;
    kdc = k*10+1;
end

%% =====================================================================================================
%% Start of Initializing Zeros ============================================
initZeros(k, kdc)

%% =====================================================================================================
%% Start Initialization ===================================================
initNLsim() % calls handleStDx at end.

%% step 3: Beginning of Huen's  (predictor-corrector) method
% Create indicies for simulation
g.k.ks = 1;
g.bus.bus_sim = g.bus.bus;

% added from v2.3 06/01/20 - thad
%g.mac.mac_trip_flags = false(g.mac.n_mac,1); - initialized during handleNewGlobals - thad 08/28/20
g.mac.mac_trip_states = 0; % seems unused...

%% ========================================================================
% Variable time step specific ==== Temporary location =====================

% initialize global st/dx vectors
handleStDx(1, 0, 0)     % init vectors name cells
handleStDx(1, [], 3)    % update g.vts.stVec to initial conditions of states
handleStDx(1, [], 1)    % update g.vts.dxVec to initial conditions of derivatives

% initlaize network solution handling
handleNetworkSln([],0)  % init network solution vector
handleNetworkSln(1, 1)  % update g.vts.netSlnVec to initial network solution

% defining ODE input and output functions
inputFcn = str2func('vtsInputFcn');
outputFcn = str2func('vtsOutputFcn');

% define ODE settings
%options = odeset('RelTol',1e-3,'AbsTol',1e-6); % MATLAB default settings
options = odeset('RelTol',1e-5,'AbsTol',1e-7, ...
    'InitialStep', 1/60/4, ...
    'MaxStep',60, ...
    'OutputFcn',outputFcn); % set 'OutputFcn' to function handle

g.vts.options = options;

% intialize counters and solution iteration log
g.vts.dataN = 1;
g.vts.iter = 0; % for keeping track of solution iterations
g.vts.tot_iter = 0;
g.vts.slns = zeros(1, size(g.sys.t,2)); % number of networks solutions used per step

% machine ref... always set to zero....?
g.sys.mach_ref = zeros(1, size(g.sys.t,2));

%% Simulation loop start
warning('*** Simulation Loop Start')
try % for catching non-convergence
    
for simTblock = 1:size(g.vts.t_block)
    
    g.vts.t_blockN = simTblock;
    g.k.ks = simTblock; % required for huen's solution method h_sol selection
    
    if ~isempty(g.vts.solver_con)
        odeName = g.vts.solver_con{g.vts.t_blockN}; % select user defined soln method
    else
        odeName = 'huens'; % default PST solver
    end
    
    % account for extra index increment from FTS to VTS - thad - 08/21/20
    if ~strcmp(odeName, 'huens') && (g.vts.t_blockN > 1)        % use VTS
        if   strcmp(g.vts.solver_con{g.vts.t_blockN-1}, 'huens')% previous was FTS
            g.vts.dataN = g.vts.dataN-1;
        end
    end
    
    % Select solution method =========================================
    if strcmp( odeName, 'huens')    % use standard PST huens method
        fprintf('*** Using Huen''s integration method for time block %d\n*** t=[%7.4f, %7.4f]\n', ...
            simTblock, g.vts.fts{simTblock}(1), g.vts.fts{simTblock}(end))
        
        % incorporate fixed time vector into system time vector
        nSteps = length(g.vts.fts{simTblock});
        g.sys.t(g.vts.dataN:g.vts.dataN+nSteps-1) = g.vts.fts{simTblock};
        
        % account for pretictor last step time check
        g.sys.t(g.vts.dataN+nSteps) = g.sys.t(g.vts.dataN+nSteps-1) + g.sys.sw_con(simTblock,7);
        
        for cur_Step = 1:nSteps
            k = g.vts.dataN;
            j = k+1;
            
            % display k and t at every first, last, and 250th step
            if ( mod(k,250)==0 ) || cur_Step == 1 || cur_Step == nSteps
                fprintf('*** k = %5d, \tt(k) = %7.4f\n',k,g.sys.t(k)) % DEBUG
                %fprintf('* current time: %s\t* start time: %s\n', datestr(now), startStr);
                % elapsed time
                %tempClock = clock;
                %eClock = tempClock - startClock;
                %fprintf('Elapsed Hours: %3d\t Minutes: %3d\t Seconds: %3.4f\n', ...
                %    eClock(4),eClock(5),eClock(6));
            end
            
            %% Time step start
            initStep(k)
            
            %% Predictor Solution =========================================
            networkSolutionVTS(k, g.sys.t(k))          
            % most recent network solution based on completely calculated states is k
            monitorSolution(k); % moved for correct AGC action -thad 08/18/20
            dynamicSolution(k)
            dcSolution(k)
            predictorIntegration(k, j, g.k.h_sol)   % g.k.h_sol updated in network solution i_simu call
            
            %% Corrector Solution =========================================
            networkSolutionVTS(j, g.sys.t(j))
            dynamicSolution(j)
            dcSolution(j)
            correctorIntegration(k, j, g.k.h_sol)
            
            %% Live plot call
            if g.sys.livePlotFlag
                livePlot(k)
            end
            
            g.vts.dataN = j;                        % increment data counter
            g.vts.tot_iter = g.vts.tot_iter  + 2;   % increment total solution counter
            g.vts.slns(g.vts.dataN) = 2;            % track step solution
        end
        % Account for next time block using VTS
        handleStDx(j, [], 3) % update g.vts.stVec to initial conditions of states
        handleStDx(k, [], 1) % update g.vts.dxVec to initial conditions of derivatives
          
    else % use user supplied variable method
        fprintf('*** Using %s integration method for time block %d\n*** t=[%7.4f, %7.4f]\n', ...
            odeName, simTblock, g.vts.t_block(simTblock, 1), g.vts.t_block(simTblock, 2))
        
        % feval used for ODE call - could be replaced with if statements.
        feval(odeName, inputFcn, g.vts.t_block(simTblock,:), g.vts.stVec , options);
        
        % Alternative example of using actual function name:
        %ode113(inputFcn, g.vts.t_block(simTblock,:), g.vts.stVec , options);
    end
    
end% end simulation loop
catch ME
    % Custom catch message
    disp('*!* Something has gone wrong and was caught!')
    fprintf('*!* Data Index: %d\t Simulation Time: %5.5f\n\n', g.vts.dataN, g.sys.t(g.vts.dataN))
    ME
    for eN = 1:size(ME.stack,1)
        revN = size(ME.stack,1)-eN+1;
        fprintf('Error in %s at line %d\n', ME.stack(revN).name, ME.stack(revN).line)
    end
end

%% ========================================================================
%% Variable step specific cleanup =========================================

% check if last time block used huens, remove last 'extra' step
if strcmp(g.vts.solver_con{end}, 'huens')
    g.vts.dataN = g.vts.dataN-1;
else
    % final network/state solution for VTS
    networkSolution(g.vts.dataN)
    dynamicSolution(g.vts.dataN)
    g.vts.tot_iter = g.vts.tot_iter + 1;
    monitorSolution(g.vts.dataN);
end

%% Display final simulation timing
disp('*** End simulation.')
et = toc;
ets = num2str(et);
g.sys.ElapsedNonLinearTime = ets;
disp(['*** Elapsed Simulation Time = ' ets 's'])

%% VTS specific output
disp(['*** Total solutions = ' int2str(g.vts.tot_iter)])
disp(['*** Total data points = ' int2str(g.vts.dataN)])

%whos('g') % to see if trimming/clearing zeros matters. 1/2

%% Trim logged values to length of g.vts.dataN
trimLogs(g.vts.dataN)

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

%% Final 'live' plot call - occurs AFTER data trimming
if g.sys.livePlotFlag
    livePlot('end')
end

%% Clean up various temporary and function input values
clear V1 V2 R X B jj phi            % used in calls to line_cur, line_pq
clear Y1 Y2 Y3 Y4 Vr1 Vr2 bo        % used in calls to i_simu, red_ybus
clear et ets                        % used to store/display elapsed time
clear nSteps fStep odeName simTblock inputFcn outputFcn options% used in VTS

clear k kdc j                       % simulation counters
clear t t_switch sw_count lswitch n_switch %used in legacy time vector creation
clear Sbase Fbase                   % used in system init
clear tol iter_max acc              % used in load flow solution

%% Remove all 'zero only' data
varNames = who()';      % all variable names in workspace
clearedVars = {};       % cell to hold names of deleted 'all zero' variables

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

%whos('g') % see if trimming/clearing zeros matters. 2/2 (it does)

%% Execute original s_simu plotting (if run as a standalone script)
standAlonePlot(scriptRunFlag)
clear scriptRunFlag

disp('*** s_simu End')
disp(' ')