function initZeros(k, kdc)
%INITZEROS Creates zero arrays for logged values
% INITZEROS Creates zero arrays for logged values based on passed in input
%
% Syntax: initZeros(k, kdc)
%
%   NOTES:
%
%   Input:
%   k - total number of time steps in the simulation
%   kdc - total number of DC time steps in the simulation
%
%   Output:
%   VOID
%
%   History:
%   Date        Time    Engineer        Description
%   07/23/20    12:54   Thad Haines     Version 1
%   08/21/20    10:54   Thad Haines     Version 1.1 - added icAdj to AGC


%% Remaining 'loose' globals
% DeltaP/omega filter variables - 21
global  dpw_con dpw_out dpw_pot dpw_pss_idx dpw_mb_idx dpw_idx n_dpw dpw_Td_idx dpw_Tz_idx
global  sdpw1 sdpw2 sdpw3 sdpw4 sdpw5 sdpw6
global  dsdpw1 dsdpw2 dsdpw3 dsdpw4 dsdpw5 dsdpw6

% pss design - 3 - Not used in Simulation? - thad 07/18/20
global ibus_con  netg_con  stab_con

%%
global g

%% create zero matrices for variables to make algorithm more efficient?
if g.sys.DEBUG
    warning('*** Initialize zero matricies')
end

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
if g.ivm.n_ivm ~= 0
    g.ivm.ivmmod_d_sig = zeros(g.ivm.n_ivm,k);
    g.ivm.ivmmod_e_sig = zeros(g.ivm.n_ivm,k);
else
    g.ivm.ivmmod_d_sig = zeros(1,k);
    g.ivm.ivmmod_e_sig = zeros(1,k);
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
        g.area.area(areaN).icAdj = zeros(1,k); % Interchange adjustment signal
        g.area.area(areaN).icS = ones(1,k); % Scheduled interchange ones for future init multiply
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

end% end initZeros