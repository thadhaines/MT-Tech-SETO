% svm_mgen.m
% 6:27 PM 18/8/97
% m.file to generate state variable models
% This m file takes the load flow and dynamic data from a data m file
% and calulates the state space matrices:
% 		A matrix in  a_mat
%               B matrices
%                 for a change in Exciter Vref in b_vr
%                 for a change in Turbine governer Pref b_pr
%                 for a change in generator mechanical torque b_pm
%                 for a change in svc reference voltage b_svc
%                 for a change in tcsc input  b_tcsc
%                 for a change in active power modulation b_lmod
%                 for a change in reactive load modulation b_rlmod
%                 for a change in real power modulation b_pwrmod_p
%                 for a change in real power modulation b_pwrmod_q
%
%               C matrices
%                 change in generator speed  c_spd
%                 change in generator electrical torque c_t on generator base
%                 change in generator electrical power c_p on generator base
%                 change in bus voltage angles c_ang
%                 change in bus voltage magnitude  c_v
%                 change in from_bus line active power c_pf1
%                 change in from_bus line reactive power c_qf1
%                 change in to_bus line active power c_pf2
%                 change in to_bus line reactive power c_qf2
%                 change in from bus current magnitude c_ilmf
%                 change in to bus current magnitude c_ilmt
%                 change in from bus real current  c_ilrf
%                 change in to bus real current c_ilrt
%                 change in from bus imaginary current c_ilif
%                 change in to bus imaginary current c_ilit
%               D matrices
%                 combination of output and input, e.g.,
%                 for power out and p_ref in --- d_ppr
% l is the eigenvalue vector of a_mat
% u is the right eigenvector matrix of a_mat
% v is the left eigenvector matrix of a_mat (vu = I)
% p is the unscaled participation matrix
% p_norm is the scaled participation matrix
% the maximum value of each column in p_norm is unity
% all scaled participations less than 0.1 are set to zero
% to find the states associated with the jth eigenvalue
% use sparse(abs(p_norm(:,j)))
% pr gives the participation factors of the eigenvalues
% associated with the rotor angle states of each generator
% use sparse(abs(pr(k,:))) to find the modes associated
% with the rotor angle of the kth generator
%
% Author: Graham Rogers
% Modified December 1998
% tcsc model added
% Modified July 1998
% deltaP/omega filter added
% Modified June 1998
% hydraulic turbine/governor added
% Modified: August 1997
% Induction Generator added
% Modified: August 1997
%           load modulation and output matrices for line flow added
% Modified: April 1997
%           HVDC added
% Version 1.0
% Date: September 1996
%
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
% clear
% clear global
tic % start timer

%% set up global variables

% DeltaP/omega filter variables
global  dpw_con dpw_out dpw_pot dpw_pss_idx dpw_mb_idx dpw_idx n_dpw dpw_Td_idx dpw_Tz_idx
global  sdpw1 sdpw2 sdpw3 sdpw4 sdpw5 sdpw6
global  dsdpw1 dsdpw2 dsdpw3 dsdpw4 dsdpw5 dsdpw6

% pss design
global  netg_con  stab_con

% ivmmod - not added to linear analysis yet - thad 09/15/20
global n_ivm mac_ivm_idx ivmmod_data ivmmod_d_sig ivmmod_e_sig

% begning of global strucutred g
global g

%% load input data from m.file
disp('linearized model development by perturbation of the non-linear model')
%set user defined SVC and TCSC models to empty
g.svc.svc_dc = [];

g.dc.dci_dc=[];
g.dc.dcr_dc=[];
% input data file
% [dfile,pathname]=uigetfile('d*.m','Select Data File');
% if pathname == 0
%    error(' you must select a valid data file')
% else
%    lfile =length(dfile);
%    % strip off .m and convert to lower case
%    dfile = lower(dfile(1:lfile-2));
%    eval(dfile);
% end
dfile = 'DataFile';
eval(dfile);

handleNewGlobals

% check for valid dynamic data file
if isempty(g.mac.mac_con)
    error(' the selected file is not a valid data file')
end

%basdat = inputdlg({'Base MVA:','Base Frequency Hz:'},'Input Base Data',1,{'100','60'});
basdat = {'100';'60'};
g.sys.sys_freq = str2double(basdat{2});
g.sys.basrad = 2*pi*g.sys.sys_freq; % default system frequency is 60 Hz
g.sys.basmva = str2double(basdat{1});
g.sys.syn_ref = 0 ;     % synchronous reference frame

% disp(' ')
% lfpf = inputdlg('do you wish to perform a post fault load flow?Y/N[N]','s');
% if isempty(lfpf{1}); lfpf{1} = 'N';end
% if lfpf{1} =='y'; lfpf{1} = 'Y'; end
% if lfpf{1} == 'Y'
%    disp('enter the changes to bus and line required to give the post fault condition')
%    disp('when you have finished, type return and press enter')
%    keyboard
% end

% solve for loadflow - loadflow parameter

if isempty(g.dc.dcsp_con)
    g.dc.n_conv = 0;
    g.dc.n_dcl = 0;
    tol = 1e-5;   % tolerance for convergence
    iter_max = 30; % maximum number of iterations
    acc = 1.0;   % acceleration factor
    [bus_sol,line,line_flw] = ...
        loadflow(g.bus.busOG, g.line.lineOG, tol,iter_max,acc,'n',2);
    g.bus.bus = bus_sol;  % solved loadflow solution needed for initialization
    g.line.line = line;
    clear bus_sol line
    %save sim_fle.mat bus line n_conv n_dcl % no need in batch runs - thad 07/17/20
else
    [bus_sol,line,line_flw,rec_par,inv_par, line_par] = lfdcs(g.bus.busOG,g.line.lineOG,g.dc.dci_dc,g.dc.dcr_dc);
    g.bus.bus = bus_sol;
    g.line.line = line;
    clear bus_sol line
    %save sim_fle.mat bus line rec_par  inv_par line_par% no need in batch runs - thad 07/17/20
end

g.exc.n_exc= 0;
g.exc.n_dc = 0;
g.exc.n_smp = 0;
g.exc.n_st3 = 0;
g.exc.n_smppi = 0; % added 07/10/20 - thad

g.pss.n_pss= 0;

n_dpw = 0;

g.tg.n_tg = 0;
g.tg.n_tgh = 0;

g.svc.n_svc = 0;
g.tcsc.n_tcsc = 0;

g.lmod.n_lmod = 0;

g.rlmod.n_rlmod = 0;

g.pwr.n_pwrmod = 0;
% note dc_indx called in load flow
mac_indx;% identifies generators
ntot = g.mac.n_mac+g.ind.n_mot+g.igen.n_ig;
ngm = g.mac.n_mac+g.ind.n_mot;
g.mac.pm_sig = zeros(g.mac.n_mac,2);
% mac_exc=0;

% check for infinite buses
if g.mac.n_ib~=0
    %remove controls associated with infinite bus generators
    %remove exciters
    if ~isempty(g.exc.exc_con)
        g.exc.n_exc = length(g.exc.exc_con(:,1));
        net_idx = zeros(g.exc.n_exc,1);
        for j = 1:g.mac.n_ib
            net_idx = net_idx | g.exc.exc_con(:,2) == g.mac.mac_con(g.mac.mac_ib_idx(j),1);
        end
        if length(net_idx)==1;
            if net_idx ==1
                g.exc.exc_con=[];
            end
        else
            perm = diag(~net_idx);
            perm = perm(~net_idx,:);
            g.exc.exc_con = perm*g.exc.exc_con;
        end
    end
    % remove pss
    if ~isempty(g.pss.pss_con)
        g.pss.n_pss = length(g.pss.pss_con(:,1));
        net_idx = zeros(g.pss.n_pss,1);
        for j = 1:g.mac.n_ib
            net_idx = net_idx | g.pss.pss_con(:,2) == g.mac.mac_con(g.mac.mac_ib_idx(j),1);
        end
        if length(net_idx)==1
            if net_idx == 1
                g.pss.pss_con = [];
            end
        else
            perm = diag(~net_idx);
            perm = perm(~net_idx,:);
            g.pss.pss_con = perm*g.pss.pss_con;
        end
    end
    % remove deltaP/omega filter
    if ~isempty(dpw_con)
        n_dpw = length(dpw_con(:,1));
        net_idx = zeros(n_dpw,1);
        for j = 1:g.mac.n_ib
            net_idx = net_idx|dpw_con(:,1) == g.mac.mac_con(g.mac.mac_ib_idx(j),1);
        end
        if length(net_idx)==1
            if net_idx == 1;dpw_con = [];end
        else
            perm = diag(~net_idx);
            perm = perm(~net_idx,:);
            dpw_con = perm*dpw_con;
        end
    end
    % remove turbine/governos
    if ~isempty(g.tg.tg_con)
        g.tg.n_tg = length(g.tg.tg_con(:,1));
        net_idx= zeros(g.tg.n_tg,1);
        for j=1:g.mac.n_ib
            net_idx =net_idx| g.tg.tg_con(:,2) == g.mac.mac_con(g.mac.mac_ib_idx(j),1);
        end
        if length(net_idx)==1
            if net_idx==1
                g.tg.tg_con = [];
            end
        else
            perm = diag(~net_idx);
            perm = perm(~net_idx,:);
            g.tg.tg_con = perm*g.tg.tg_con;
        end
    end
end

if ~isempty(g.exc.exc_con)
    exc_indx;%identifies exciters
    mac_exc = g.mac.mac_int(g.exc.exc_con(:,2));
else
    g.exc.n_exc=0;
    mac_exc=0; % added -thad 07/10/20
end

mac_pss=0;
if ~isempty(g.pss.pss_con)
    pss_indx;%identifies power system stabilizers
    mac_pss= g.mac.mac_int(g.pss.pss_con(:,2));
else
    g.pss.n_pss=0;
end
mac_dpw=0;
if ~isempty(dpw_con)
    dpwf_indx;%identifies deltaP/omega filters
    mac_dpw= g.mac.mac_int(dpw_con(:,1));
else
    n_dpw=0;
end

if ~isempty(g.tg.tg_con)
    tg_indx;%identifies turbine/governor
    mac_tg = g.mac.mac_int(g.tg.tg_con(g.tg.tg_idx,2));
    mac_tgh = g.mac.mac_int(g.tg.tg_con(g.tg.tgh_idx,2));
else
    g.tg.n_tg =0;
    g.tg.n_tgh = 0;
    mac_tg = 0; % added 07/10/20
    mac_tgh = 0;
end

if ~isempty(g.svc.svc_con)~=0
    g.svc.svc_dc=[];
    svc_indx();
else
    g.svc.n_svc = 0;
end

g.tcsc.tcsc_dc=[];
g.tcsc.n_tcscud=0;

if ~isempty(g.tcsc.tcsc_con)
    tcsc_indx();
else
    g.tcsc.n_tcsc = 0;
end

if ~isempty(g.lmod.lmod_con)
    lm_indx; % identifies load modulation buses
    % line flow monitoring buses? (Chow, 02/28/2016)
else
    g.lmod.n_lmod = 0;
end

if ~isempty(g.rlmod.rlmod_con)~=0
    rlm_indx; % identifies load modulation buses
else
    g.rlmod.n_rlmod = 0;
end

if ~isempty(g.pwr.pwrmod_con)
    pwrmod_indx(g.bus.bus); % identifies power modulation buses % corrected to call pwrmod_indx, not pwrm_indx
else
    g.pwr.n_pwrmod = 0;
end

%initialize induction motor
if g.ind.n_mot~=0
    g.ind.vdp = zeros(g.ind.n_mot,2);
    g.ind.vqp = zeros(g.ind.n_mot,2);
    g.ind.slip = zeros(g.ind.n_mot,2);
    g.ind.dvdp = zeros(g.ind.n_mot,2);
    g.ind.dvqp = zeros(g.ind.n_mot,2);
    g.ind.dslip = zeros(g.ind.n_mot,2);
end
g.bus.bus = mac_ind(0,1,g.bus.bus,0);

%initialize induction generator
if g.igen.n_ig~=0
    g.igen.vdpig = zeros(g.igen.n_ig,2);
    g.igen.vqpig = zeros(g.igen.n_ig,2);
    g.igen.slig = zeros(g.igen.n_ig,2);
    g.igen.dvdpig = zeros(g.igen.n_ig,2);
    g.igen.dvqpig = zeros(g.igen.n_ig,2);
    g.igen.dslig = zeros(g.igen.n_ig,2);
    g.igen.tmig = zeros(g.igen.n_ig,2);
end
g.bus.bus = mac_igen(0,1,g.bus.bus,0);

%initialize svc
if g.svc.n_svc ~=0
    g.svc.B_cv = zeros(g.svc.n_svc,2);
    g.svc.dB_cv = zeros(g.svc.n_svc,2);
    g.svc.B_con = zeros(g.svc.n_svc,2);
    g.svc.dB_con = zeros(g.svc.n_svc,2);
    if g.svc.n_dcud~=0
        error('user defined svc damping control not allowed in small signal simulation')
    else
        g.svc.svc_dsig = zeros(g.svc.n_svc,2);
    end
end
g.bus.bus = svc(0,1,g.bus.bus,0);

if g.dc.n_conv~=0
    % pick up HVDC initial variables from load flow
    g.dc.Vdc(g.dc.r_idx,1) = rec_par(:,2);
    g.dc.Vdc(g.dc.i_idx,1) = inv_par(:,2);
    g.dc.i_dc(g.dc.r_idx,1) = line_par;
    g.dc.i_dc(g.dc.i_idx,1) = line_par;
    g.dc.i_dcr(:,1) = g.dc.i_dc(g.dc.r_idx,1);
    g.dc.i_dci(:,1) = g.dc.i_dc(g.dc.i_idx,1);
    g.dc.alpha(:,1) = rec_par(:,1)*pi/180;
    g.dc.gamma(:,1) = inv_par(:,1)*pi/180;
    g.dc.dcr_dsig = zeros(g.dc.n_dcl,2);
    g.dc.dci_dsig = zeros(g.dc.n_dcl,2);
    g.dc.dc_sig = zeros(g.dc.n_conv,2);
end
dc_cont(0,1,1,g.bus.bus,0); % initialize the dc controls - sets up data for red_ybus
% this has to be done before red_ybus is used since the induction motor,svc and
% dc link initialization alters the bus matrix
v = ones(length(g.bus.bus(:,1)),2);
g.bus.bus_v = v;
g.bus.theta = zeros(length(g.bus.bus(:,1)),2);
disp(' ')
disp('Performing linearization')
% set line parameters
if ~isempty(g.lmon.lmon_con)
    R = g.line.line(g.lmon.lmon_con,3);
    X = g.line.line(g.lmon.lmon_con,4);
    B = g.line.line(g.lmon.lmon_con,5);
    g.dc.tap = g.line.line(g.lmon.lmon_con,6); % seems odd.... -thad 07/14/20
    phi = g.line.line(g.lmon.lmon_con,7);
end
% step 1: construct reduced Y matrix
[Y_gprf,Y_gncprf,Y_ncgprf,Y_ncprf,V_rgprf,V_rncprf,boprf] = red_ybus(g.bus.bus,g.line.line);
g.bus.bus_intprf = g.bus.bus_int;% store the internal bus numbers for the pre_fault system
nbus = length(g.bus.bus(:,1));
if isempty(g.ncl.load_con)
    g.ncl.nload = 0;
else
    g.ncl.nload = length(g.ncl.load_con(:,1));
end
state = zeros(g.mac.n_mac,1);
gen_state = state;
TR_state = state;
TB_state = state;
TA_state = state;
Efd_state = state;
R_f_state = state;
pss1_state = state;
pss2_state = state;
pss3_state = state;
dpw_state = state;
tg_state = state;
state = zeros(g.mac.n_mac+g.ind.n_mot+g.igen.n_ig+g.svc.n_svc+g.tcsc.n_tcsc ...
    +g.lmod.n_lmod + g.rlmod.n_rlmod+2*g.pwr.n_pwrmod+g.dc.n_dcl,1);
max_state = 6*g.mac.n_mac + 5*g.exc.n_exc+ 3*g.pss.n_pss+ 6*n_dpw ...
    + 5*g.tg.n_tg+ 5*g.tg.n_tgh+ 3*g.ind.n_mot+ 3*g.igen.n_ig+ ...
    2*g.svc.n_svc+g.tcsc.n_tcsc+ g.lmod.n_lmod  +g.rlmod.n_rlmod+2*g.pwr.n_pwrmod+5*g.dc.n_dcl;
%25 states per generator,3 per motor, 3 per ind. generator,
% 2 per SVC,1 per tcsc, 1 per lmod,1 per rlmod, 2 per pwrmod, 5 per dc line
g.bus.theta(:,1) = g.bus.bus(:,3)*pi/180;
v(:,1) = g.bus.bus(:,2).*exp(1j*g.bus.theta(:,1));
if g.dc.n_conv ~= 0
    % convert dc LT to Equ HT bus
    Pr = g.bus.bus(g.dc.rec_ac_bus,6);
    Pi = g.bus.bus(g.dc.inv_ac_bus,6);
    Qr = g.bus.bus(g.dc.rec_ac_bus,7);
    Qi = g.bus.bus(g.dc.inv_ac_bus,7);
    VLT= v(g.dc.ac_bus,1);
    i_acr = (Pr-1j*Qr)./conj(VLT(g.dc.r_idx));
    i_aci = (Pi - 1j*Qi)./conj(VLT(g.dc.i_idx));
    v(g.dc.rec_ac_bus,1) = VLT(g.dc.r_idx) + 1j*g.dc.dcc_pot(:,2).*i_acr;
    v(g.dc.inv_ac_bus,1) = VLT(g.dc.i_idx) + 1j*g.dc.dcc_pot(:,4).*i_aci;
    g.bus.theta(g.dc.ac_bus,1) = angle(v(g.dc.ac_bus,1));
end
g.bus.bus_v(:,1) = v(:,1);
v(:,2) = v(:,1);
g.bus.bus_v(:,2)=v(:,1);
g.bus.theta(:,2) = g.bus.theta(:,1);
% find total number of states
ns_file
NumStates = sum(state);
g.exc.exc_sig = zeros(g.mac.n_mac,2);

if n_dpw~=0; dpw_out = zeros(n_dpw,2);else dpw_out = zeros(1,2);end
if g.tg.n_tg ~=0||g.tg.n_tgh ~= 0
    g.tg.tg_sig = zeros(g.tg.n_tg+g.tg.n_tgh,2);
else
    g.tg.tg_sig = zeros(1,2);
end
if g.svc.n_svc ~=0
    g.svc.svc_sig = zeros(g.svc.n_svc,2);
else
    g.svc.svc_sig = zeros(1,2);
end
if g.tcsc.n_tcsc ~=0
    g.tcsc.tcsc_sig = zeros(g.tcsc.n_tcsc,2);
else
    g.tcsc.tcsc_sig = zeros(1,2);
end
if g.lmod.n_lmod ~= 0
    g.lmod.lmod_sig = zeros(g.lmod.n_lmod,2);
else
    g.lmod.lmod_sig = zeros(1,2);
end
if g.rlmod.n_rlmod ~= 0
    g.rlmod.rlmod_sig = zeros(g.rlmod.n_rlmod,2);
else
    g.rlmod.rlmod_sig = zeros(1,2);
end
if g.pwr.n_pwrmod ~= 0
    g.pwr.pwrmod_p_sig = zeros(g.pwr.n_pwrmod,2);
    g.pwr.pwrmod_q_sig = zeros(g.pwr.n_pwrmod,2);
else
    g.pwr.pwrmod_p_sig = zeros(1,2);
    g.pwr.pwrmod_q_sig = zeros(1,2);
end

if g.dc.n_conv ~= 0
    g.dc.dc_sig = zeros(g.dc.n_conv,2);
else
    g.dc.dc_sig = zeros(2,2);
end
% set initial state and rate matrices to zero
nMacZero = zeros(g.mac.n_mac,2);
g.mac.psi_re = nMacZero;
g.mac.psi_im = nMacZero;
psi = nMacZero;

g.mac.eterm = nMacZero;
g.mac.pelect = nMacZero;
g.mac.qelect = nMacZero;
g.mac.mac_ang = nMacZero;
g.mac.mac_spd = nMacZero;
g.mac.edprime = nMacZero;
g.mac.eqprime = nMacZero;
g.mac.psikd = nMacZero;
g.mac.psikq = nMacZero;
g.mac.dmac_ang = nMacZero;
g.mac.dmac_spd = nMacZero;
g.mac.dedprime = nMacZero;
g.mac.deqprime = nMacZero;
g.mac.dpsikd = nMacZero;
g.mac.dpsikq = nMacZero;

if g.exc.n_exc~=0
    g.exc.V_TR = zeros(g.exc.n_exc,2);
    g.exc.V_As = zeros(g.exc.n_exc,2);
    g.exc.V_A = zeros(g.exc.n_exc,2);
    g.exc.V_R =zeros(g.exc.n_exc,2);
    g.exc.Efd = zeros(g.exc.n_exc,2);
    g.exc.R_f = zeros(g.exc.n_exc,2);
    g.exc.dV_TR = zeros(g.exc.n_exc,2);
    g.exc.dV_As = zeros(g.exc.n_exc,2);
    g.exc.dV_R =zeros(g.exc.n_exc,2);
    g.exc.dEfd = zeros(g.exc.n_exc,2);
    g.exc.dR_f = zeros(g.exc.n_exc,2);
    
    g.pss.pss_out = zeros(g.exc.n_exc,2);
end

if g.pss.n_pss~=0
    g.pss.pss1 = zeros(g.pss.n_pss,2);
    g.pss.pss2 = zeros(g.pss.n_pss,2);
    g.pss.pss3 = zeros(g.pss.n_pss,2);
    g.pss.dpss1 = zeros(g.pss.n_pss,2);
    g.pss.dpss2 =zeros(g.pss.n_pss,2);
    g.pss.dpss3 =zeros(g.pss.n_pss,2);
end
if n_dpw~=0
    sdpw1 = zeros(n_dpw,2);
    sdpw2 = zeros(n_dpw,2);
    sdpw3 = zeros(n_dpw,2);
    sdpw4 = zeros(n_dpw,2);
    sdpw5 = zeros(n_dpw,2);
    sdpw6 = zeros(n_dpw,2);
    dsdpw1 = zeros(n_dpw,2);
    dsdpw2 = zeros(n_dpw,2);
    dsdpw3 = zeros(n_dpw,2);
    dsdpw4 = zeros(n_dpw,2);
    dsdpw5 = zeros(n_dpw,2);
    dsdpw6 = zeros(n_dpw,2);
end
if g.tg.n_tg~=0 || g.tg.n_tgh~=0
    tgZeros = zeros(g.tg.n_tg+g.tg.n_tgh,2);
    g.tg.tg1 = tgZeros;
    g.tg.tg2 = tgZeros;
    g.tg.tg3 = tgZeros;
    g.tg.tg4 = tgZeros;
    g.tg.tg5 = tgZeros;
    g.tg.dtg1 = tgZeros;
    g.tg.dtg2 = tgZeros;
    g.tg.dtg3 = tgZeros;
    g.tg.dtg4 = tgZeros;
    g.tg.dtg5 = tgZeros;
    clear tgZeros
end
if g.lmod.n_lmod~=0
    g.lmod.lmod_st = zeros(g.lmod.n_lmod,2);
    g.lmod.dlmod_st = zeros(g.lmod.n_lmod,2);
end
if g.rlmod.n_rlmod~=0
    g.rlmod.rlmod_st = zeros(g.rlmod.n_rlmod,2);
    g.rlmod.drlmod_st = zeros(g.rlmod.n_rlmod,2);
end
if g.pwr.n_pwrmod~=0
    g.pwr.pwrmod_p_st = zeros(g.pwr.n_pwrmod,2);
    g.pwr.dpwrmod_p_st = zeros(g.pwr.n_pwrmod,2);
    g.pwr.pwrmod_q_st = zeros(g.pwr.n_pwrmod,2);
    g.pwr.dpwrmod_q_st = zeros(g.pwr.n_pwrmod,2);
end

%HVDC links
if g.dc.n_conv~= 0
    g.dc.dv_conr = zeros(g.dc.n_dcl,2);
    g.dc.dv_coni = zeros(g.dc.n_dcl,2);
    g.dc.di_dcr = zeros(g.dc.n_dcl,2);
    g.dc.di_dci = zeros(g.dc.n_dcl,2);
    g.dc.dv_dcc = zeros(g.dc.n_dcl,2);
end
% set dimensions for A matrix and permutation matrix
a_mat = zeros(NumStates);
p_mat = sparse(zeros(NumStates,max_state));
c_spd = zeros(length(g.mac.not_ib_idx),NumStates);
c_p = zeros(length(g.mac.not_ib_idx),NumStates);
c_t = zeros(length(g.mac.not_ib_idx),NumStates);

%determine p_mat: converts the vector of length max_states to
%a column of a_mat or b

p_m_file

% step 2: initialization
flag = 0;
%machines
mac_em(0,1,g.bus.bus,flag);
mac_tra(0,1,g.bus.bus,flag);
mac_sub(0,1,g.bus.bus,flag);
mac_ib(0,1,g.bus.bus,flag);
%calculate initial electrical torque
psi = g.mac.psi_re(:,1)+1j*g.mac.psi_im(:,1);
if g.ind.n_mot~=0&&g.igen.n_ig==0
    vmp = g.ind.vdp(:,1) + 1j*g.ind.vqp(:,1);
    int_volt=[psi; vmp]; % internal voltages of generators and motors
elseif g.ind.n_mot==0&&g.igen.n_ig~=0
    vmpig = g.igen.vdpig(:,1) + 1j*g.igen.vqpig(:,1);
    int_volt = [psi; vmpig]; % int volt of synch and ind generators
elseif g.ind.n_mot~=0&&g.igen.n_ig~=0
    vmp = g.ind.vdp(:,1) + 1j*g.ind.vqp(:,1);
    vmpig = g.igen.vdpig(:,1) + 1j*g.igen.vqpig(:,1);
    int_volt = [psi; vmp; vmpig];
else
    int_volt = psi;
end
cur(:,1) = Y_gprf*int_volt; % network solution currents into generators
b_v(boprf(g.ncl.nload+1:nbus),1) = V_rgprf*int_volt;   % bus voltage reconstruction
if g.ncl.nload~=0
    vnc = nc_load(g.bus.bus,flag,Y_ncprf,Y_ncgprf);%vnc is a dummy variable
    cur(:,1) = cur(:,1) + Y_gncprf*v(g.bus.bus_intprf(g.ncl.load_con(:,1)),1);% modify currents for nc loads
end
g.mac.cur_re(1:g.mac.n_mac,1) = real(cur(1:g.mac.n_mac,1));
g.mac.cur_im(1:g.mac.n_mac,1) = imag(cur(1:g.mac.n_mac,1));
cur_mag(1:g.mac.n_mac,1) = abs(cur(1:g.mac.n_mac,1)).*g.mac.mac_pot(:,1);
if g.ind.n_mot~=0
    g.ind.idmot(:,1) = -real(cur(g.mac.n_mac+1:ngm,1));%induction motor currents
    g.ind.iqmot(:,1) = -imag(cur(g.mac.n_mac+1:ngm,1));%current out of network
end
if g.igen.n_ig~=0
    g.igen.idig(:,1) = -real(cur(ngm+1:ntot,1));%induction generator currents
    g.igen.iqig(:,1) = -imag(cur(ngm+1:ntot,1));%current out of network
end

if g.dc.n_conv ~=0
    % calculate dc voltage and current
    V0(g.dc.r_idx,1) = abs(v(g.dc.rec_ac_bus,1)).*g.dc.dcc_pot(:,7);
    V0(g.dc.i_idx,1) = abs(v(g.dc.inv_ac_bus,1)).*g.dc.dcc_pot(:,8);
    g.dc.Vdc(g.dc.r_idx,1) = V0(g.dc.r_idx,1).*cos(g.dc.alpha(:,1)) - g.dc.i_dcr(:,1).*g.dc.dcc_pot(:,3);
    g.dc.Vdc(g.dc.i_idx,1) = V0(g.dc.i_idx,1).*cos(g.dc.gamma(:,1)) - g.dc.i_dci(:,1).*g.dc.dcc_pot(:,5);
    g.dc.Vdc_ref = g.dc.Vdc(g.dc.i_idx,1);
    g.dc.i_dc(g.dc.r_idx,1) = g.dc.i_dcr(:,1);
    g.dc.i_dc(g.dc.i_idx,1) = g.dc.i_dci(:,1);
end

telect(:,1) = g.mac.pelect(:,1).*g.mac.mac_pot(:,1)...
    + cur_mag(:,1).*cur_mag(:,1).*g.mac.mac_con(:,5);
% DeltaP/omega filter
dpwf(0,1,flag);
%pss
pss(0,1,flag);
%exciters
smpexc(0,1,flag);
smppi(0,1,flag);
exc_dc12(0,1,flag);
exc_st3(0,1,flag);
% turbine governors
tg(0,1,flag);
tg_hydro(0,1,g.bus.bus,flag);

%initialize tcsc
if g.tcsc.n_tcsc ~=0
    g.tcsc.B_tcsc = zeros(g.tcsc.n_tcsc,2);
    g.tcsc.dB_tcsc = zeros(g.tcsc.n_tcsc,2);
    if g.tcsc.n_tcscud~=0
        error('user defined tcsc damping control not allowed in small signal simulation')
    else
        g.tcsc.tcsc_dsig = zeros(g.tcsc.n_tcsc,2);
    end
end
tcsc(0,1,0);

if ~isempty(g.lmod.lmod_con)
    disp('load modulation')
    lmod(0,1,flag);
end
if ~isempty(g.rlmod.rlmod_con)
    disp('reactive load modulation')
    rlmod(0,1,flag);
end
if ~isempty(g.pwr.pwrmod_con)
    disp('power modulation')
    pwrmod_p(0,1,g.bus.bus,flag);
    pwrmod_q(0,1,g.bus.bus,flag);
end

% initialize non-linear loads
if ~isempty(g.ncl.load_con)
    disp('non-linear loads')
    vnc = nc_load(g.bus.bus,flag,Y_ncprf,Y_ncgprf);
else
    g.ncl.nload = 0;
end
% hvdc lines
dc_line(0,1,1,g.bus.bus,flag);

g.sys.mach_ref(1) = 0;
g.sys.mach_ref(2) = 0;
g.sys.sys_freq(1) = 1;
g.sys.sys_freq(2) = 1;

%set states
%generators
g.mac.mac_ang(:,2) = g.mac.mac_ang(:,1);
g.mac.mac_spd(:,2) = g.mac.mac_spd(:,1);
g.mac.eqprime(:,2) = g.mac.eqprime(:,1);
g.mac.psikd(:,2) = g.mac.psikd(:,1);
g.mac.edprime(:,2) = g.mac.edprime(:,1);
g.mac.psikq(:,2)= g.mac.psikq(:,1);

%exciters
if g.exc.n_exc~=0
    g.exc.V_TR(:,2)=g.exc.V_TR(:,1);
    g.exc.V_As(:,2) = g.exc.V_As(:,1);
    g.exc.V_A(:,2) = g.exc.V_A(:,1);
    g.exc.V_R(:,2)=g.exc.V_R(:,1);
    g.exc.Efd(:,2)=g.exc.Efd(:,1);
    g.exc.R_f(:,2)=g.exc.R_f(:,1);
end

%pss
if g.pss.n_pss~=0
    g.pss.pss1(:,2)=g.pss.pss1(:,1);
    g.pss.pss2(:,2)=g.pss.pss2(:,1);
    g.pss.pss3(:,2)=g.pss.pss3(:,1);
end
% DeltaP/omega filter
if n_dpw~=0
    sdpw1(:,2)=sdpw1(:,1);
    sdpw2(:,2)=sdpw2(:,1);
    sdpw3(:,2)=sdpw3(:,1);
    sdpw4(:,2)=sdpw4(:,1);
    sdpw5(:,2)=sdpw5(:,1);
    sdpw6(:,2)=sdpw6(:,1);
end
%turbine governor
if g.tg.n_tg~=0 || g.tg.n_tgh~=0
    g.tg.tg1(:,2) = g.tg.tg1(:,1);
    g.tg.tg2(:,2) = g.tg.tg2(:,1);
    g.tg.tg3(:,2) = g.tg.tg3(:,1);
    g.tg.tg4(:,2) = g.tg.tg4(:,1);
    g.tg.tg5(:,2) = g.tg.tg5(:,1);
end

telect(:,2) =telect(:,1); % unused? -thad 07/09/20
if g.ind.n_mot~=0
    g.ind.vdp(:,2) = g.ind.vdp(:,1);
    g.ind.vqp(:,2) = g.ind.vqp(:,1);
    g.ind.slip(:,2) = g.ind.slip(:,1);
end
if g.igen.n_ig~=0
    g.igen.vdpig(:,2) = g.igen.vdpig(:,1);
    g.igen.vqpig(:,2) = g.igen.vqpig(:,1);
    g.igen.slig(:,2) = g.igen.slig(:,1);
    g.igen.tmig(:,2) = g.igen.tmig(:,1);
end
if g.svc.n_svc ~= 0
    g.svc.B_cv(:,2) = g.svc.B_cv(:,1);
    g.svc.B_con(:,2) = g.svc.B_con(:,1);
end
if g.tcsc.n_tcsc ~= 0
    g.tcsc.B_tcsc(:,2) = g.tcsc.B_tcsc(:,1);
end
if g.lmod.n_lmod ~=0
    g.lmod.lmod_st(:,2) = g.lmod.lmod_st(:,1);
end
if g.rlmod.n_rlmod ~=0
    g.rlmod.rlmod_st(:,2) = g.rlmod.rlmod_st(:,1);
end
if g.pwr.n_pwrmod ~=0
    g.pwr.pwrmod_p_st(:,2) = g.pwr.pwrmod_p_st(:,1);
    g.pwr.pwrmod_q_st(:,2) = g.pwr.pwrmod_q_st(:,1);
end
if g.dc.n_conv~=0
    g.dc.v_conr(:,2) = g.dc.v_conr(:,1);
    g.dc.v_coni(:,2) = g.dc.v_coni(:,1);
    g.dc.i_dcr(:,2) = g.dc.i_dcr(:,1);
    g.dc.i_dci(:,2) = g.dc.i_dci(:,1);
    g.dc.v_dcc(:,2) = g.dc.v_dcc(:,1);
end
% set interconnection variables in perturbation stage to defaults
% this accounts for any generators which do not have the
% corresponding controls
g.mac.eterm(:,2) = g.mac.eterm(:,1);
g.mac.pmech(:,2) = g.mac.pmech(:,1);
g.mac.vex(:,2) = g.mac.vex(:,1);
g.exc.exc_sig(:,2) = g.exc.exc_sig(:,1);
g.tg.tg_sig(:,2) = g.tg.tg_sig(:,1);
g.svc.svc_sig(:,2) = g.svc.svc_sig(:,1);
g.tcsc.tcsc_sig(:,2) = g.tcsc.tcsc_sig(:,1);
g.lmod.lmod_sig(:,2) = g.lmod.lmod_sig(:,1);
g.rlmod.rlmod_sig(:,2) = g.rlmod.rlmod_sig(:,1);
if g.pwr.n_pwrmod ~=0
    g.pwr.pwrmod_p_sig(:,1) = g.pwr.pwrmod_p_st(:,1);
    g.pwr.pwrmod_q_sig(:,1) = g.pwr.pwrmod_q_st(:,1);
end
g.pwr.pwrmod_p_sig(:,2) = g.pwr.pwrmod_p_sig(:,1);
g.pwr.pwrmod_q_sig(:,2) = g.pwr.pwrmod_q_sig(:,1);

if g.dc.n_conv~=0
    g.dc.Vdc(:,2) = g.dc.Vdc(:,1);
    g.dc.i_dc(:,2) = g.dc.i_dc(:,1);
    g.dc.dc_sig(:,2) = g.dc.dc_sig(:,1);
    g.dc.cur_ord(:,2) = g.dc.cur_ord(:,1);
    g.dc.alpha(:,2) = g.dc.alpha(:,1);
    g.dc.gamma(:,2) = g.dc.gamma(:,1);
end
%perform perturbation of state variables

p_cont

% setup matrix giving state numbers for generators
mac_state = zeros(sum(state(1:g.mac.n_mac)),3);
for k = 1:g.mac.n_mac
    if state(k)~=0
        if k == 1
            j = 1;
        else
            j = 1+sum(state(1:k-1));
        end
        jj = sum(state(1:k));
        mac_state(j:jj,1) = (j:jj)';
        mac_state(j:jj,2) = st_name(k,1:state(k))';
        mac_state(j:jj,3) = k*ones(state(k),1);
    end
end

ang_idx = find(mac_state(:,2)==1);
b_pm = zeros(NumStates,g.mac.n_mac-g.mac.n_ib);

b_pm(ang_idx+1,:)= diag(0.5./g.mac.mac_con(g.mac.not_ib_idx,16));
% Form transformation matrix to get rid of zero eigenvalue
% Use generator 1 as reference
% check for infinite buses
if isempty(g.mac.ibus_con)
    ref_gen = 'N';%questdlg('Set gen 1 as reference');
    if strcmp(ref_gen,'Yes')
        p_ang = eye(NumStates);
        p_ang(ang_idx,1) = -ones(length(ang_idx),1);
        p_ang(1,1) = 1;
        p_angi = inv(p_ang);
        %transform state matrix
        a_mat = p_ang*a_mat*p_angi;
        %transform the c matrices
        c_v = c_v*p_angi;
        c_ang = c_ang*p_angi;
        c_spd = c_spd*p_angi;
        c_pm = c_pm*p_angi;
        c_t = c_t*p_angi;
        c_p = c_p*p_angi;
        if ~isempty(g.lmon.lmon_con)
            c_pf1 = c_pf1*p_angi;
            c_qf1 = c_qf1*p_angi;
            c_pf2 = c_pf2*p_angi;
            c_qf2 = c_qf2*p_angi;
            c_ilmf = c_ilmf*p_angi;
            c_ilmt = c_ilmt*p_angi;
            c_ilrf = c_ilrf*p_angi;
            c_ilrt = c_ilrt*p_angi;
            c_ilif = c_ilif*p_angi;
            c_ilit = c_ilit*p_angi;
        end
        if g.dc.n_conv~=0;
            c_dcir=c_dcir*p_angi;
            c_dcii=c_dcii*p_angi;
            c_Vdcr=c_dcVr*p_angi;
            c_Vdci=c_Vdci*p_angi;
        end
        %transform the b matrices
        b_pm = p_ang*b_pm;
        if g.tg.n_tg~=0
            b_pr = p_ang*b_pr;
        end
        if g.exc.n_exc~=0
            b_vr = p_ang*b_vr;
        end
        if g.svc.n_svc~=0
            b_svc = p_ang*b_svc;
        end
        if g.tcsc.n_tcsc~=0
            b_tcsc = p_ang*b_tcsc;
        end
        if g.lmod.n_lmod~=0
            b_lmod = p_ang*b_lmod;
        end
        if g.rlmod.n_rlmod~=0
            b_rlmod = p_ang*b_rlmod;
        end
        if g.pwr.n_pwrmod~=0
            b_pwrmod_p = p_ang*b_pwrmod_p;
        end
        if g.pwr.n_pwrmod~=0
            b_pwrmod_q = p_ang*b_pwrmod_q;
        end
        if g.exc.n_dc~=0
            b_dcr=p_ang*b_dcr;
            b_dci=p_ang*b_dci;
        end
    end
end

disp('calculating eigenvalues and eigenvectors')
%eigenvectors and eigenvalues of a_mat
[u l] = eig(a_mat); % u is the right eigenvector

% sort the eigenvalues
[l l_idx] =sort( diag(l));


%reorder the eigenvector matrix
u = u(:,l_idx);

for j = 1:NumStates
    if imag(l(j))~=0
        %scale the complex eigenvectors so that the maximum element is 1+j0
        [maxu,mu_idx] = max(abs(u(:,j)));
        u(:,j) = u(:,j)/u(mu_idx,j);
    end
end
v = inv(u); % left eigenvectors
% find the participation factors
p=zeros(NumStates);p_norm=p;
for j = 1:NumStates
    p(:,j) = (conj(v(j,:)))'.*u(:,j);% p are the unnormalized participation vectors
    [p_max,p_max_idx] = max((p(:,j)));
    p_norm(:,j) = p(:,j)/p(p_max_idx,j);% p_norm has biggest element = 1
    p_big = abs(p_norm(:,j))>0.1;%big sorts out normalized participation > 1
    p_norm(:,j) = p_big.*p_norm(:,j);% p_norm now contains only values of p_norm > 0.1
end

% find states associated with the generator angles
pr = p_norm(ang_idx,:);
% frequency and damping ratio
freq = abs(imag(l))/2/pi;
zero_eig = find(abs(l)<=1e-4);
if ~isempty(zero_eig)
    damp(zero_eig,1)= ones(length(zero_eig),1);
end
nz_eig = find(abs(l)>1e-4);
damp(nz_eig,1) = -real(l(nz_eig))./abs(l(nz_eig));

% figure
% hold on
% box on
% stab_idx =find(damp>=0.05);
% plot(damp(stab_idx),freq(stab_idx),'k+')
% fmax = ceil(max(freq));
% plot([0.05 0.05],[0 fmax],'r')
% title(['Calculated Modes ' dfile])
% xlabel('damping ratio')
% ylabel('frequency Hz')
% us_idx = find(damp<0);
% plot(damp(us_idx),freq(us_idx),'r+')
% ud_idx = find(damp>0&damp<0.05);
% plot(damp(ud_idx),freq(ud_idx),'g+')

% full sim timing
et = toc;
ets = num2str(et);
g.sys.ElapsedLinearTime = ets;
disp(['elapsed time = ' ets 's'])
disp('*** End simulation.')
disp(' ')

%% tidy work space
% clear global
clear B Efd_state Pi Pr Qi Qr R R_f_state
clear TA_state TB_state TR_state V0 V1 V2 VLT V_rgprf V_rncprf
clear X Y_gncprf Y_gprf Y_ncgprf Y_ncprf
clear ans b_v boprf bus_intprf bvnc
clear c_state chk_dc chk_smp cur cur_mag
clear d_vector dc_start dci_dc dcmod_input
clear dfile  dpw_count dpw_state exc_count exc_number exc_state
clear f   flag from_idx  gen_state  gh
clear i_aci  i_acr int_volt inv_par
clear j  j_state jay jj k k_cex k_col k_colg k_ctg k_dc k_exc k_exc_idx
clear k_hvdc k_idx k_lmod k_nib_idx k_rlmod k_row k_smp k_sub k_tg k_tgh k_tra  kgs
clear l_idx l_if1 l_if2 l_it1 l_it2 lf lfile
clear line_flw line_par lmod_input lmod_start mac_dpw mac_exc mac_pss mac_tg mac_tgh
clear max_state maxu mu_idx n_hvdc1 n_hvdc_states n_ig_states n_lmod1 n_lmod_states
clear n_mot_states n_rlmod1 n_rlmod_states n_svc_states n_tcsc_states nbus ngm ngt
clear no_mac nominal not_TA not_TB not_TE not_TF not_TR not_ib
clear ntdc  ntf  ntl ntot ntrl nts nz_eig ntpwr_p ntpwr_q
clear p_ang p_angi p_big p_mat p_max p_max_idx p_ratio pathname pert
clear phi pr pr_input  psi pss1_state pss2_state pss3_state  pss_count  pss_state
clear rec_par ref_gen rlmod_input rlmod_start  s11 s12  s21  s22
clear s_TA s_TB s_TE s_TR s_idx sel st_name
clear state  state_hvdc state_lmod state_rlmod
clear telect tg_count tg_number tg_state to_idx vnc vr_input zero_eig
