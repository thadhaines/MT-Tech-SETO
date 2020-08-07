function predictorIntegration(k, j, h_sol)
%PREDICTORINTEGRATION Performs x(j) = x(k) + h_sol*dx(k)
% PREDICTORINTEGRATION Performs  x(j) = x(k) + h_sol*dx(k)
%
% Syntax: predictorIntegration(k, j, h_sol)
%
%   NOTES:  
%
%   Input:
%   k - data index for 'n'
%   j - data index for 'n+1'
%   h_sol - time between k and j
%
%   Output:
%   VOID
%
%   History:
%   Date        Time    Engineer        Description
%   07/23/20    13:10   Thad Haines     Version 1
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
                g.pwr.pwrmod_p_sigst{index}(:,j) = g.pwr.pwrmod_p_sigst{index}(:,k)+h_sol*g.pwr.dpwrmod_p_sigst{index}(:,k);
                g.pwr.pwrmod_q_sigst{index}(:,j) = g.pwr.pwrmod_q_sigst{index}(:,k)+h_sol*g.pwr.dpwrmod_q_sigst{index}(:,k);
            end
        end
        
        %% ivmmod
        if n_ivm~=0
            for index=1:n_ivm
                ivmmod_d_sigst{index}(:,j) = ivmmod_d_sigst{index}(:,k)+h_sol*divmmod_d_sigst{index}(:,k);
                ivmmod_e_sigst{index}(:,j) = ivmmod_e_sigst{index}(:,k)+h_sol*divmmod_e_sigst{index}(:,k);
            end
        end
        
        %% agc predictor integration
        if g.agc.n_agc ~=0
            for ndx = 1:g.agc.n_agc
                g.agc.agc(ndx).sace(j) = g.agc.agc(ndx).sace(k) + h_sol*g.agc.agc(ndx).d_sace(k);
                
                % integrate lowpass filter outs...
                for gndx=1:g.agc.agc(ndx).n_ctrlGen
                    g.agc.agc(ndx).ctrlGen(gndx).x(j) = g.agc.agc(ndx).ctrlGen(gndx).x(k)...
                        + h_sol * g.agc.agc(ndx).ctrlGen(gndx).dx(k)  ;
                end
            end
        end