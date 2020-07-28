function handleStDx(k, slnVec, flag)
%HANDLESTDX Performs required state and derivative handling for ODE solvers
% HANDLESTDX Performs various state and derivative handling functions to
% allow for MATLAB varialbe timestep integration techniques. 
% Uses dynamic field names. 
%
% Syntax: handleStDx(k, slnVec, flag)
%
%   NOTES:  Requires state and derivative values are in the same g.(x) field.
%           Not all flags require same input/output
%
%   Input:
%   k - data index
%   flag - choose between operations
%           0 - initialize state and derivative cell array, count states
%           1 - update g.vts.dxVec with col k of derivative fields
%           2 - write slnVec vector of values to associated states at index k+1
%           3 - return vector of values to associated states at index k
%           4 - DEBUG verify write (kind of untested/unused)
%-
%   snlVec - used to populated states with new values
%
%   Output:
%   VOID
%
%   History:
%   Date        Time    Engineer        Description
%   07/24/20    10:45   Thad Haines     Version 1
%   07/27/20    09:57   Thad Haines     Version 1.0.1 - integrated into VTS sim

%% Remaining 'loose' globals - Not required for now -thad 07/24/20
% % ivm variables - 5
% global n_ivm mac_ivm_idx ivmmod_data ivmmod_d_sig ivmmod_e_sig
%
% % DeltaP/omega filter variables - 21
% global  dpw_con dpw_out dpw_pot dpw_pss_idx dpw_mb_idx dpw_idx n_dpw dpw_Td_idx dpw_Tz_idx
% global  sdpw1 sdpw2 sdpw3 sdpw4 sdpw5 sdpw6
% global  dsdpw1 dsdpw2 dsdpw3 dsdpw4 dsdpw5 dsdpw6
%
% % pss design - 3 - Not used in Simulation? - thad 07/18/20
% global ibus_con  netg_con  stab_con

%%
global g

if flag == 0
    
    %% Initialize Field, State, Derivative, Number cell
    % no  input, no output
    % action: create cell array of states and derivaties,
    %   count number of states, init global dxVec
    
    fsdn = { ...
        %  Field 1, State Field 2,  Derivative Field 2, nStates
        % machine
        'mac', 'mac_ang', 'dmac_ang', 0 ;
        'mac', 'mac_spd', 'dmac_spd', 0 ;
        'mac', 'edprime', 'dedprime', 0 ;
        'mac', 'eqprime', 'deqprime', 0 ;
        'mac', 'psikd', 'dpsikd', 0 ;
        'mac', 'psikq', 'dpsikq', 0 ;
        % exciter
        'exc', 'Efd', 'dEfd', 0;
        'exc', 'V_R', 'dV_R', 0;
        'exc', 'V_As', 'dV_As', 0;
        'exc', 'R_f', 'dR_f', 0;
        'exc', 'V_TR', 'dV_TR', 0;
        % pss
        'pss', 'pss1', 'dpss1', 0;
        'pss', 'pss2', 'dpss2', 0;
        'pss', 'pss3', 'dpss3', 0;
        % governor
        'tg', 'tg1', 'dtg1', 0;
        'tg', 'tg2', 'dtg2', 0;
        'tg', 'tg3', 'dtg3', 0;
        'tg', 'tg4', 'dtg4', 0;
        'tg', 'tg5', 'dtg5', 0;
        };
    
    % induction motor
    if g.ind.n_mot ~= 0
        fsdn = [fsdn; {'ind', 'vdp','dvdp',0}];
        fsdn = [fsdn; {'ind', 'vqp','dvqp',0}];
        fsdn = [fsdn; {'ind', 'slip','dslip',0}];
    end
    
    % induction generator
    if g.igen.n_ig ~=0
        fsdn = [fsdn; {'igen', 'vdpig','dvdpig',0}];
        fsdn = [fsdn; {'igen', 'vqpig','dvqpig',0}];
        fsdn = [fsdn; {'igen', 'slig','dslig',0}];
    end
    
    % svc
    if g.svc.n_svc ~= 0
        fsdn = [fsdn; {'svc', 'B_cv','dB_cv',0}];
        fsdn = [fsdn; {'svc', 'B_con','dB_con',0}];
        fsdn = [fsdn; {'svc', 'xsvc_dc','dxsvc_dc',0}];
    end
    
    %tcsc
    if g.tcsc.n_tcsc ~= 0
        fsdn = [fsdn; {'tcsc', 'B_tcsc','dB_tcsc',0}];
        fsdn = [fsdn; {'tcsc', 'xtcsc_dc','dxtcsc_dc',0}];
    end
    
    if g.lmod.n_lmod~=0
        fsdn = [fsdn; {'lmod', 'lmod_st','dlmod_st',0}];
    end
    if g.rlmod.n_rlmod~=0
        fsdn = [fsdn; {'rlmod', 'rlmod_st','drlmod_st',0}];
    end
    
    %fsdn % debug output
    
    % Count number of allocated states, update fsdn
    nStates = 0;
    for ndx=1:size(fsdn,1)
        f1 = fsdn{ndx,1};
        f2 = fsdn{ndx, 2};
        fsdn{ndx, 4} = size( g.(f1).(f2), 1);
        nStates = nStates + fsdn{ndx, 4};
    end
    
    % store required globals g
    g.vts.fsdn = fsdn; 
    g.vts.n_states = nStates;
    g.vts.dxVec = zeros(nStates,1);
    g.vts.stVec = zeros(nStates,1);
    
    %fsdn % debug output
end

if flag == 1
    %% update g.vts.dxVec with col k of derivative fields
    % input: index k to collect derivatives from
    % output: VOID
    % action: updated dxVec with derivatives from index k
    
    %fprintf('DEBUG: Creating single vector of all g.x.dxstate(:,k)...\n')
    
    startN = 1; % used to split up vector according to number of states
    
    for ndx=1:size(g.vts.fsdn,1)
        
        % If the number of allocated states is larger than zero
        if g.vts.fsdn{ndx, 4} ~= 0
            % collect derivative loction
            f1 = g.vts.fsdn{ndx,1};
            f2 = g.vts.fsdn{ndx, 3};
            
            % place derivatives into dxVec using dynamic field names
            g.vts.dxVec(startN:startN+g.vts.fsdn{ndx,4}-1) = g.(f1).(f2)(1:g.vts.fsdn{ndx,4}, k);
            
            startN = startN + g.vts.fsdn{ndx,4}; % increment starting dxVec index
        else
            % probably not executed...
            disp(['skpping... ', g.vts.fsdn(ndx,1), g.vts.fsdn(ndx,3) ] );
        end
    end
    
    %dxVec % debug output
end

if flag == 2
    %% Write slnVec vector of values to states at index k+1
    % input =   slnVec - internally generated by the ODE sovler
    %           k - index of derivative states. Writes next states (integration results) to k+1
    % output = none
    
    %fprintf('DEBUG: Performing g.x.state(:,j) = g.x.dxstate(:,k)...\n')
    
    j = k+1;
    startN = 1;
    for ndx=1:size(g.vts.fsdn,1)
        % If the number of allocated states is larger than zero
        if g.vts.fsdn{ndx, 4} ~= 0
            % collect field and state name
            f1 = g.vts.fsdn{ndx,1};
            f2 = g.vts.fsdn{ndx, 2};
            for i=1:g.vts.fsdn{ndx, 4}
                g.(f1).(f2)(i,j) = slnVec(startN+i-1); % write to state col j
            end
            startN = startN+g.vts.fsdn{ndx, 4};
        else
            % probably not executed...
            disp(['skpping... ', g.vts.fsdn(ndx,1), g.vts.fsdn(ndx,2) ] );
        end
    end
    
end

if flag == 3
    %% update state vector of states at time index k
    % input k
    % output, none
    % action - updated g.vts.stVec (used for intialization of ODE sovler)
    
    startN = 1; % used to split up vector according to number of states
    
    for ndx=1:size(g.vts.fsdn,1)        
        % If the number of allocated states is larger than zero
        if g.vts.fsdn{ndx, 4} ~= 0
            % collect derivative loction
            f1 = g.vts.fsdn{ndx, 1};
            f2 = g.vts.fsdn{ndx, 2};
            
            % place state into stVec using dynamic field names
            g.vts.stVec(startN:startN+g.vts.fsdn{ndx,4}-1) = g.(f1).(f2)(1:g.vts.fsdn{ndx,4}, k);
            
            startN = startN + g.vts.fsdn{ndx,4}; % increment starting stVec index
        else
            % probably not executed...
            disp(['skpping... ', g.vts.fsdn(ndx,1), g.vts.fsdn(ndx,3) ] );
        end
    end
    
end

if flag == 4
    %% DEBUG verify write
    % input ...
    % include at all?
    for ndx=1:size(g.vts.fsdn,1)
        if g.vts.fsdn{ndx, 4} ~= 0
            % collect field and state name
            f1 = g.vts.fsdn{ndx,1};
            f2 = g.vts.fsdn{ndx, 2};
            f3 = g.vts.fsdn{ndx, 3};
            disp([f1, '.',f2,'(:,j)',' to ', f1,'.',f3 '(:,k)'] );
            for i=1:g.vts.fsdn{ndx, 4}
                v1 = g.(f1).(f2)(i,j);
                v2 = g.(f1).(f3)(i,k);
                fprintf('Check: %10.7f\t%10.7f\n',v1,v2)
            end
        else
            disp(['skpping... ', g.vts.fsdn(ndx,1), g.vts.fsdn(ndx,2) ] );
        end
    end
end

%% unique/unused state and derivatives not currently handled in function
%{
% Copied from v2.3 - 06/01/20 - thad
g.pwr.pwrmod_p_st(:,j) = g.pwr.pwrmod_p_st(:,k)+h_sol*(g.pwr.dpwrmod_p_st(:,j) + g.pwr.dpwrmod_p_st(:,k))/2;
g.pwr.pwrmod_q_st(:,j) = g.pwr.pwrmod_q_st(:,k)+h_sol*(g.pwr.dpwrmod_q_st(:,j) + g.pwr.dpwrmod_q_st(:,k))/2;
if g.pwr.n_pwrmod~=0
    for index=1:g.pwr.n_pwrmod
        % make global? -thad 07/06/20
        pwrmod_p_sigst{index}(:,j) = pwrmod_p_sigst{index}(:,k)+h_sol*(dpwrmod_p_sigst{index}(:,j) + dpwrmod_p_sigst{index}(:,k))/2;
        pwrmod_q_sigst{index}(:,j) = pwrmod_q_sigst{index}(:,k)+h_sol*(dpwrmod_q_sigst{index}(:,j) + dpwrmod_q_sigst{index}(:,k))/2;
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

if n_dpw ~= 0
    % only calculate if dpw filter is used
    sdpw1(:,j) = sdpw1(:,k) +h_sol*(dsdpw1(:,k)+dsdpw1(:,j))/2.;
    sdpw2(:,j) = sdpw2(:,k) +h_sol*(dsdpw2(:,k)+dsdpw2(:,j))/2.;
    sdpw3(:,j) = sdpw3(:,k) +h_sol*(dsdpw3(:,k)+dsdpw3(:,j))/2.;
    sdpw4(:,j) = sdpw4(:,k) +h_sol*(dsdpw4(:,k)+dsdpw4(:,j))/2.;
    sdpw5(:,j) = sdpw5(:,k) +h_sol*(dsdpw5(:,k)+dsdpw5(:,j))/2.;
    sdpw6(:,j) = sdpw6(:,k) +h_sol*(dsdpw6(:,k)+dsdpw6(:,j))/2.;
end
%}