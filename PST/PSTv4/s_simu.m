%S_SIMU runs a non-linear simulation of a given data file.
% S_SIMU simulates power system transients using Power System Toolbox
% This m-file takes the dynamic and load flow data and calculates the
% response of the power system to a fault (which is specified in a
% switching file) or other power system perturbance. See one of the supplied
% examples for file format and/or replacing technique.
%
%   NOTES:  To run in stand-alone mode, clear all variables before
%           running script from PST main directory (i.e. executre 'clear all').
%           Accomodates for batch mode runs automatically if a non-global
%           variable is detected.
%
%   See user manual for details and flow chart of code operation
%
%   Input:
%   N/A
%
%   Output:
%   g - structured gloabl variable containing trimmed and cleaned data
%
%   History:
%{
Date        Time    Engineer        Description
02/xx/97    XX:XX   Graham Rogers	Version 1.0
08/xx/97    xx:xx   Graham Rogers   Version 1.1 - add modulation for load,
                                    exciter, and governor.
08/xx/97    XX:XX   Graham Rogers   Version 1.2 - add induction generator
(c) Copyright: Joe Chow/ Cherry Tree Scientific Software 1991 to 1997
                                    - All rights reserved
06/xx/98    xx:xx   Graham Rogers   Version 1.3 - Add hydraulic turbine/governor
07/xx/98    XX:XX   Graham Rogers   Version 1.4 - add deltaP/omega filter
12/xx/98    xx:xx   Graham Rogers   Version 1.5 - add svc user defined
                                    damping controls.
01/xx/99    XX:XX   Graham Rogers   Version 1.5 - add tcsc model and tcsc user
                                    defined damping controls.
06/14/99    09:59   Graham Rogers   Version 1.6 - add user defined hvdc control
                                    modify dc so that it is integrated with a
                                    sub-multiple time constant
09/xx/99    XX:XX   Graham Rogers   Version 1.7 - Add simple exciter with
                                    pi avr - smppi
xx/xx/15    xx:xx   Dan Trudnowski  Version 1.75 - Added pwrmod code
07/19/17    xx:xx   Rui             Version 1.8 - Add code for multiple HVDC
xx/xx/19    xx:xx   Dan Trudnowski  Version 1.9 - Added ivmmod code
07/15/20    14:27   Thad Haines     Version 2.0 - Revised format of globals
                                    and internal function documentation
07/16/20    11:19   Thad Haines     V 2.1 - Added cleanZeros to end of script
                                    to clean global g
07/18/20    10:52   Thad Haines     V 2.2 - Added Line monitor, area, and
                                    sytem average frequency calculations.
07/21/20    16:20   Thad Haines     V 2.3 - Added AGC
07/23/20    11:24   Thad Haines     V 2.3.1 - Begining work on
                                    functionalization of solution.
07/27/20    09:54   Thad Haines     V 2.3.2 - Begining work on VTS simulation
07/28/20    16:20   Thad Haines     V 2.4 - Initial working VTS simulation
07/29/20    15:20   Thad Haines     V 2.4.1 - jay -> 1j
07/30/20    08:04   Thad Haines     V 2.5 - Re-integrated interactive script
                                    running/plotting
07/30/20    20:34   Thad Haines     V 2.6 - Added selectable solution methods
08/11/20    11:36   Thad Haines     V 2.7 - added ivm to global
08/13/20    09:37   Thad Haines     V 2.8 - Incorporated new license
08/16/20    20:28   Thad Haines     V 2.8.1 - Increase logged data only if VTS
08/18/20    10:51   Thad Haines     V 2.8.2 - Moved network solution in Huen's
                                    method to work with AGC
08/21/20    12:53   Thad Haines     V 2.8.3 - Handled FTS->VTS unique time
                                    vector issue
08/25/20    11:02   Thad Haines     V 2.8.4 - Added try catch to handle
                                    non-convergence
08/28/20    14:05   Thad Haines     V 2.8.5 - removed mac_trip_flags init,
                                    changed number of steps to print data
09/08/20    04:41   Thad Haines     V 2.8.6 - added option to not solve a
                                    power flow - to stand alone mode
09/12/20    10:13   Thad Haines     V 2.8.7 - code clean up, addition of
                                    g.sys.DEBUG flag
10/01/20    09:43   Thad Haines     4.0.0 recheck - dub
10/13/20    11:25   Thad Haines     4.0.0 release
10/19/20    13:52   Thad Haines     4.0.1 - fix of area_indx
%}

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
disp('***    PST v4.0.1    ***')
disp('***')
disp('*** s_simu Start')
disp('*** Declaring Global Variables')

%% Remaining 'loose' globals - UNUSED (probably...)
%% DeltaP/omega filter variables - model use not documented
global dpw_con dpw_out dpw_pot dpw_pss_idx dpw_mb_idx
global dpw_idx n_dpw dpw_Td_idx dpw_Tz_idx
global sdpw1 sdpw2 sdpw3 sdpw4 sdpw5 sdpw6
global dsdpw1 dsdpw2 dsdpw3 dsdpw4 dsdpw5 dsdpw6

%% pss design -  Not used  - thad 07/18/20
global ibus_con  netg_con  stab_con

%% global structured array
global g
% set debug flag for addition display
g.sys.DEBUG = 0;

%% Check for Octave, automatically load compatibility script
% Assumes license of Octave will be 'GNU ...'
dataCheck = license;
if all(dataCheck(1:3)=='GNU')
    fprintf(['*** Octave detected,',...
        ' loading compatiblity commands and packages...\n'])
    octaveComp
else
    fprintf('*** MATLAB detected\n')
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
    copyfile([pathname, fname],'DataFile.m'); % copy system data file
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
if ~exist('Fbase','var')
    fprintf('*** Fbase not defined - assuming 60 Hz base.\n')
    g.sys.sys_freq = 60;
    g.sys.Fbase = 60;
elseif isnumeric(Fbase)
    fprintf('*** Fbase found - Frequency base is set to %3.3f Hz\n', Fbase)
    g.sys.sys_freq = Fbase;
    g.sys.Fbase = Fbase;
end
g.sys.basrad = 2*pi*g.sys.sys_freq;

% Handle variable input base MVA
if ~exist('Sbase','var')
    fprintf('*** Sbase not defined - assuming 100 MVA base.\n')
    g.sys.basmva = 100;
elseif isnumeric(Sbase)
    fprintf('*** Sbase found - Power base is set to %3.3f MVA\n', Sbase)
    g.sys.basmva = Sbase;
end

clear Fbase Sbase

%% unused/unimplemented damping controls?
% intentionally removed/ignored? -thad 07/15/20
g.svc.svc_dc = [];
g.tcsc.tcsc_dc = [];
g.tcsc.n_tcscud = 0;
g.dc.dcr_dc = [];
g.dc.dci_dc = [];
g.sys.syn_ref = 0 ;     % synchronous reference frame
g.mac.ibus_con = [];    % ignore infinite buses in transient simulation

%% solve for loadflow - loadflow parameter
solveLoadFlow = 1;
disp('***')

if scriptRunFlag == 1
    % ask if load flow should be solved
    prompt={'Solve Load Flow? [ y / n ]' };
    name='Should a Load Flow be computed?';
    numlines=1;
    defaultanswer={'y'};
    
    answer=inputdlg(prompt,name,numlines,defaultanswer);
    
    if ~strcmp(answer, 'y')
        solveLoadFlow = 0;
    end
    clear answer prompt name numlines defaultanswer answer
end

if solveLoadFlow
    if g.sys.DEBUG
        warning('*** Solve Initial Loadflow')
    else
        disp('*** Solve Initial Loadflow')
    end
    if isempty(g.dc.dcsp_con)
        % AC power flow
        g.dc.n_conv = 0;
        g.dc.n_dcl = 0;
        g.dc.ndcr_ud = 0;
        g.dc.ndci_ud = 0;
        tol = 1e-9;   % tolerance for convergence
        iter_max = 30; % maximum number of iterations
        acc = 1.0;   % acceleration factor
        [bus_sol,line,~] = ...
            loadflow(g.bus.busOG,g.line.lineOG,tol,iter_max,acc,'n',2);
        g.bus.bus = bus_sol;  % solved loadflow solution req for initialization
        g.line.line = line;
        clear bus_sol line tol iter_max acc
    else
        % Has HVDC, use DC load flow
        [bus_sol,line,~,rec_par, inv_par, line_par] = ...
            lfdcs(g.bus.busOG,g.line.lineOG,g.dc.dci_dc,g.dc.dcr_dc);
        g.bus.bus = bus_sol;
        g.line.line = line;
        g.dc.rec_par = rec_par;
        g.dc.inv_par = inv_par;
        g.dc.line_par = line_par;
        clear bus_sol line rec_par inv_par line_par
    end
else
    % Don't perform a power flow.
    g.dc.n_conv = 0;
    g.dc.n_dcl = 0;
    g.dc.ndcr_ud = 0;
    g.bus.bus = g.bus.busOG;  % solved loadflow solution req for initialization
    g.line.line = g.line.lineOG;
end
clear solveLoadFlow
disp('***')

%% Timer start
% (after all user input options)
tic

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

%% construct simulation switching sequence as defined in sw_con (OG method)
% included for reference only
if g.sys.DEBUG
    warning('*** Initialize time and switching variables')
end
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
    g.k.k_inc(sw_count) = fix((g.sys.sw_con(sw_count+1,1) ...
        -g.sys.sw_con(sw_count,1))/g.k.h(sw_count));%nearest lower integer
    
    if g.k.k_inc(sw_count) == 0
        g.k.k_inc(sw_count) = 1;
    end% minimum 1
    
    % adjust time step so integer number of steps in block
    g.k.h(sw_count) = (g.sys.sw_con(sw_count+1,1) ...
        -g.sys.sw_con(sw_count,1))/g.k.k_inc(sw_count);%step length
    g.k.h_dc(sw_count) = g.k.h(sw_count)/10; % h_dc 10 times faster
    g.k.k_incdc(sw_count) = 10*g.k.k_inc(sw_count);
    
    t_switch(sw_count+1) = t_switch(sw_count) ...
        + g.k.k_inc(sw_count)*g.k.h(sw_count);
    
    % create time vector block from start time to end time - 1 timestep
    t(k:k-1+g.k.k_inc(sw_count)) = ...
        t_switch(sw_count):g.k.h(sw_count):t_switch(sw_count+1)...
        - g.k.h(sw_count);
    
    if ~isempty(g.dc.dcl_con)
        g.dc.t_dc(kdc:kdc-1+g.k.k_incdc(sw_count)) = ...
            t_switch(sw_count):g.k.h_dc(sw_count):t_switch(sw_count+1) ...
            - g.k.h_dc(sw_count);
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

% k is the total number of steps for simulation (+1 for final predictor step)
k = sum(g.k.k_inc)+1;
t(k) = g.sys.sw_con(n_switch,1); % final time into time vector

% NOTE: the time vector created above is used only to create zeros...
g.sys.t = zeros(1,size(t,2));
g.sys.t_OLD = t; % for reference only

%% creation of VTS time blocks =================================================
initTblocks()

% If variable step, increase amount of zeros to log
if ~all(strcmp(g.vts.solver_con, 'huens'))
    k = max(size(g.sys.t_OLD))*4; % create 4x FTS amount of zeros
    g.sys.t = zeros(1,k);
    kdc = k*10;
else
    k = max(size(g.sys.t))+1;
    kdc = k*10+1;
end

%% Initialize Zeros ============================================================
initZeros(k, kdc)

%% Initialize Simulation =======================================================
initNLsim()

%% step 3: Beginning of Huen's  (predictor-corrector) method
% Create indicies for simulation
g.k.ks = 1;
g.bus.bus_sim = g.bus.bus;

% added from v2.3 06/01/20 - thad
g.mac.mac_trip_states = 0; % seems unused...

%% Variable time step specific ==== 'Temporary' location =======================

% initialize global st/dx vectors
handleStDx(1, 0, 0)     % init vectors name cells
handleStDx(1, [], 3)    % update g.vts.stVec to init cond of states
handleStDx(1, [], 1)    % update g.vts.dxVec to init cond of derivatives

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
g.vts.slns = zeros(1, size(g.sys.t,2)); % for num of network solns per step

% machine ref... always set to zero....?
g.sys.mach_ref = zeros(1, size(g.sys.t,2));

disp('*** Initialization Complete')
disp('***')

%% Simulation Loop
if g.sys.DEBUG
    warning('*** Simulation Loop Start')
else
    disp('*** Simulation Loop Start')
end

try % used to catch non-convergence or other issues/bugs
    for simTblock = 1:size(g.vts.t_block)
        
        g.vts.t_blockN = simTblock; % set required index
        g.k.ks = simTblock; % required for huen's method h_sol selection
        
        odeName = g.vts.solver_con{g.vts.t_blockN}; % select solution method
        
        % if current time block is VTS and not the first
        if ~strcmp(odeName, 'huens') && (g.vts.t_blockN > 1)
            % if previous time block was FTS
            if   strcmp(g.vts.solver_con{g.vts.t_blockN-1}, 'huens')
                % remove extra index increment for FTS to VTS time blocks
                g.vts.dataN = g.vts.dataN-1;
            end
        end
        
        %% Select solution method ==============================================
        if strcmp( odeName, 'huens')    % use standard PST Huen's method
            huensMethod(simTblock)
        else
            % use VTS
            fprintf(['*~~ Using %s integration method ',...
                'for time block %d\n*** t_block = [%7.4f, %7.4f ]\n'], ...
                odeName, simTblock, g.vts.t_block(simTblock, 1), ...
                g.vts.t_block(simTblock, 2))
            
            % feval used for ODE call - could be replaced with if statements.
            feval(odeName, inputFcn, g.vts.t_block(simTblock,:), ...
                g.vts.stVec , options);
            
            % NOTE: While feval is not recommended for repeated use,
            % the funtion is called only once per time block.
            % Example of using ODE function (ode113) name:
            % ode113(inputFcn, g.vts.t_block(simTblock,:), g.vts.stVec, options)
        end
        
    end% end simulation loop
catch ME
    % Custom catch message display
    disp('*!* Something has gone wrong and was caught!')
    fprintf('*!* Data Index: %d\t Simulation Time: %5.5f\n\n', ...
        g.vts.dataN, g.sys.t(g.vts.dataN))
    disp()
    ME % display caught error
    disp()
    for eN = 1:size(ME.stack,1)
        revN = size(ME.stack,1)-eN+1;
        fprintf('*!* Error in %s at line %d\n', ...
            ME.stack(revN).name, ME.stack(revN).line)
    end
end

%% Post Simulation Variable Step Specific Cleanup ==============================

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
disp('*** Simulation Loop End')
disp('***')
et = toc;
ets = num2str(et);
g.sys.ElapsedNonLinearTime = ets;
disp(['*** Elapsed Simulation Time = ' ets 's'])

%% VTS specific output
disp(['*** Total Solutions = ' int2str(g.vts.tot_iter)])
disp(['*** Total Data Points = ' int2str(g.vts.dataN)])

%whos('g') % to see if trimming/clearing zeros matters. 1/2

%% Trim logged values to length of g.vts.dataN
trimLogs(g.vts.dataN)

%% Clean up logged DC variables to length of DC simulated time.
% Note: could be included in trimLogs - thad 09/08/20
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
clear nSteps fStep odeName simTblock inputFcn outputFcn options % used in VTS

clear k kdc j                       % simulation counters
clear t t_switch sw_count lswitch n_switch %used in legacy time vector creation
clear Sbase Fbase                   % used in system init
clear tol iter_max acc              % used in load flow solution

%% Remove all 'zero only' data
varNames = who()';      % all variable names in workspace
clearedVars = {};       % cell to hold names of deleted 'all zero' variables

for vName = varNames
    try
        % check if all zeros
        zeroTest = eval(sprintf('all(%s(:)==0)', vName{1}));
        if zeroTest
            eval(sprintf('clear %s',vName{1}) ); % clear variable
            clearedVars{end+1} = vName{1}; % add name to cell for reference
        end
    catch %ME
        % catches stuctures
        testStr = ['isstruct(',vName{1},')'];
        if eval(testStr)
            fprintf('*** Clearing zeros from structure %s...\n', vName{1});
            eval(sprintf('[%s, clearedVars] = cleanZeros(%s, clearedVars);', ...
                vName{1},vName{1} ))
        end
        clear testStr
    end
    
end
g.sys.clearedVars = clearedVars; % attach cleard vars to global g
clear varNames vName zeroTest clearedVars % associated with clearing zeros.

%whos('g') % see if trimming/clearing zeros matters. 2/2 (it does)

%% Execute original s_simu plotting (if run as a stand alone script)
standAlonePlot(scriptRunFlag)
clear scriptRunFlag

disp('*** s_simu End')
disp(' ')