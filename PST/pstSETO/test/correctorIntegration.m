function correctorIntegration(k, j, h_sol)
%CORRECTORINTEGRATION Performs x(j) = x(k) + h(dx(j) + dx(k))/2
% CORRECTORINTEGRATION Performs x(j) = x(k) + h(dx(j) + dx(k))/2
%
% Syntax: correctorIntegration(k, j, h)
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


        %% following statements are corrector steps (RK2 computation)
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
                g.pwr.pwrmod_p_sigst{index}(:,j) = g.pwr.pwrmod_p_sigst{index}(:,k)+h_sol*(g.pwr.dpwrmod_p_sigst{index}(:,j) + g.pwr.dpwrmod_p_sigst{index}(:,k))/2;
                g.pwr.pwrmod_q_sigst{index}(:,j) = g.pwr.pwrmod_q_sigst{index}(:,k)+h_sol*(g.pwr.dpwrmod_q_sigst{index}(:,j) + g.pwr.dpwrmod_q_sigst{index}(:,k))/2;
            end
        end
        
        if n_ivm~=0
            for index=1:n_ivm
                % make global? -thad 07/06/20
                ivmmod_d_sigst{index}(:,j) = ivmmod_d_sigst{index}(:,k)+h_sol*(divmmod_d_sigst{index}(:,j) + divmmod_d_sigst{index}(:,k))/2;
                ivmmod_e_sigst{index}(:,j) = ivmmod_e_sigst{index}(:,k)+h_sol*(divmmod_e_sigst{index}(:,j) + divmmod_e_sigst{index}(:,k))/2;
            end
        end
        
        %% agc corrector integration
        if g.agc.n_agc ~=0
            for ndx = 1:g.agc.n_agc
                g.agc.agc(ndx).sace(j) = g.agc.agc(ndx).sace(k) + h_sol*(g.agc.agc(ndx).d_sace(j)+g.agc.agc(ndx).d_sace(k))/2;
                
                % integrate lowpass filter outs...
                for gndx=1:g.agc.agc(ndx).n_ctrlGen
                    g.agc.agc(ndx).ctrlGen(gndx).x(j) = g.agc.agc(ndx).ctrlGen(gndx).x(k)...
                        + h_sol *(g.agc.agc(ndx).ctrlGen(gndx).dx(j)+  g.agc.agc(ndx).ctrlGen(gndx).dx(k))/2  ;
                end
            end
        end
        