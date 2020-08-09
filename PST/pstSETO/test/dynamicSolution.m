function dynamicSolution(k)
%DYNAMICSOLUTION Performs the dynamic solution for index k
% DYNAMICSOLUTION Performs the dynamic solution for index k
%
% Syntax: dynamicSolution(k)
%
%   NOTES:  Non-global IVM variables will NOT work until added the
%           global g.
%
%   Input:
%   k - data index
%
%   Output:
%   VOID
%
%   History:
%   Date        Time    Engineer        Description
%   07/23/20    11:31   Thad Haines     Version 1
%   08/05/20    13:15   Thad Haines     Version 1.0.1 - added pwrmod signals to global

%% Remaining 'loose' globals
% ivm variables - 5
global n_ivm mac_ivm_idx ivmmod_data ivmmod_d_sig ivmmod_e_sig

% DeltaP/omega filter variables - 21
global  dpw_con dpw_out dpw_pot dpw_pss_idx dpw_mb_idx dpw_idx n_dpw dpw_Td_idx dpw_Tz_idx
global  sdpw1 sdpw2 sdpw3 sdpw4 sdpw5 sdpw6
global  dsdpw1 dsdpw2 dsdpw3 dsdpw4 dsdpw5 dsdpw6

% pss design - 3 - Not used in Simulation? - thad 07/18/20
global ibus_con  netg_con  stab_con

%%
global g

%% step 3b: compute dynamics and integrate
flag = 2;
%g.sys.sys_freq(k) = 1.0; % why?... 5/21/20 -thad
% initialized as all ones

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
        Pst{index} = g.pwr.pwrmod_p_sigst{index}(:,k);
        Qst{index} = g.pwr.pwrmod_q_sigst{index}(:,k);
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
        if size(dp{index},1)~=size(g.pwr.dpwrmod_p_sigst{index},1)
            error('Dimension error in pwrmod_dyn');
        end
        if size(dq{index},1)~=size(g.pwr.dpwrmod_q_sigst{index},1)
            error('Dimension error in pwrmod_dyn');
        end
        g.pwr.dpwrmod_p_sigst{index}(:,k) = dp{index};
        g.pwr.dpwrmod_q_sigst{index}(:,k) = dq{index};
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
end % end dynamicSolution