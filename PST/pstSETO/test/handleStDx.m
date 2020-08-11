function handleStDx(k, slnVec, flag)
%HANDLESTDX Performs required state and derivative handling for ODE solvers
% HANDLESTDX Performs various state and derivative handling functions to
% allow for MATLAB varialbe timestep integration techniques.
% Uses dynamic field names.
%
% Syntax: handleStDx(k, slnVec, flag)
%
%   NOTES:  Requires state and derivative values are in the same g.(x) field.
%           Not all flags require same input.
%
%   Input:
%   k - data index
%   flag - choose between operations
%           0 - initialize state and derivative cell array, count states
%           1 - update g.vts.dxVec with col k of derivative fields
%           2 - write slnVec vector of values to associated states at index k
%           3 - update g.vts.stVec with col k of state fields
%           4 - DEBUG verify write (kind of untested/unused)
%   snlVec - Input used to populate states with new values
%
%   Output:
%   VOID
%
%   History:
%   Date        Time    Engineer        Description
%   07/24/20    10:45   Thad Haines     Version 1
%   07/27/20    09:57   Thad Haines     Version 1.0.1 - integrated into VTS sim
%   08/03/20    15:31   Thad Haines     Version 1.0.2 - added AGC
%   08/05/20    13:15   Thad Haines     Version 1.0.3 - added pwrmod
%   08/11/20    11:18   Thad Haines     Version 1.0.4 - added ivmmod


%% Remaining 'loose' globals - Not required for now -thad 07/24/20

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
    
    % real load modulation
    if g.lmod.n_lmod~=0
        fsdn = [fsdn; {'lmod', 'lmod_st','dlmod_st',0}];
    end
    
    % reactive load modulation
    if g.rlmod.n_rlmod~=0
        fsdn = [fsdn; {'rlmod', 'rlmod_st','drlmod_st',0}];
    end
    
    % AGC
    % Each AGC has a sace and d_sace state and derivative
    % and a variable number of controlled generator x and dx
    if g.agc.n_agc ~=0
        fsdn = [fsdn; {'agc', {'sace','x'}, {'d_sace','dx'}, 0}];
    end
    
    % PWRMOD
    % each pwrmod has p, q state as well as
    % a variable number of p and q 'signal states'
    if g.pwr.n_pwrmod~=0
        fsdn = [fsdn; {'pwr','pwrmod_p_st', 'dpwrmod_p_st', 0}];
        fsdn = [fsdn; {'pwr','pwrmod_q_st', 'dpwrmod_q_st', 0}];
        % the below values are in cells
        fsdn = [fsdn; {'pwr','pwrmod_p_sigst', 'dpwrmod_p_sigst', 0}];
        fsdn = [fsdn; {'pwr','pwrmod_q_sigst', 'dpwrmod_q_sigst', 0}];
    end
    
    % IVMMOD
    % each pwrmod has d and e signal states in a cell
    if g.ivm.n_ivm ~= 0
        fsdn = [fsdn; {'ivm','ivmmod_d_sigst', 'divmmod_d_sigst', 0}];
        fsdn = [fsdn; {'ivm','ivmmod_e_sigst', 'divmmod_e_sigst', 0}];
    end
    
    %fsdn % debug output
    
    % Count number of allocated states, update fsdn
    nStates = 0;
    for ndx=1:size(fsdn,1)
        f1 = fsdn{ndx, 1};
        f2 = fsdn{ndx, 2};
        
        if strcmp(f1, 'agc') % agc handling
            agcStateVec = zeros(g.agc.n_agc,1); % each row will be number of states associated with each agc
            agcStateTot = 0;
            % for each agc
            for n = 1: g.agc.n_agc
                % collect number of agc
                agcStateVec(n) = 1 + g.agc.agc(n).n_ctrlGen; % SACE, ctrl_gen LPF
                agcStateTot = agcStateTot + agcStateVec(n);
            end
            fsdn{ndx, 4} = agcStateVec;
            nStates = nStates + agcStateTot;
            
        else
            fsdn{ndx, 4} = size( g.(f1).(f2), 1);
            nStates = nStates + fsdn{ndx, 4};
        end
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
        if all(g.vts.fsdn{ndx, 4} ~= 0)
            % collect derivative loction
            f1 = g.vts.fsdn{ndx, 1};
            f2 = g.vts.fsdn{ndx, 3};
            
            if strcmp(f1, 'agc')
                for n = 1:g.agc.n_agc
                    % collect d_sace for each agc
                    g.vts.dxVec(startN) = g.agc.agc(n).d_sace(k);
                    startN = startN+1;
                    for cg = 1:g.agc.agc(n).n_ctrlGen
                        % collect low pass dx for each controlled gen
                        g.vts.dxVec(startN) = g.agc.agc(n).ctrlGen(cg).dx(k);
                        startN = startN +1;
                    end
                end
                
            elseif strcmp(f2, 'dpwrmod_p_sigst') || strcmp(f2, 'dpwrmod_q_sigst')  || ...  % pwrmod handling of variable cell states
                   strcmp(f2, 'divmmod_d_sigst') || strcmp(f2, 'divmmod_e_sigst') % ivmmod handling of variable cell states
                for n = 1:g.vts.fsdn{ndx,4}
                    % place derivatives into dxVec using dynamic field names
                    g.vts.dxVec(startN) = g.(f1).(f2){n}(k);
                    startN = startN + 1; % increment starting dxVec index
                end
                
            else% 'normal' operation
                % place derivatives into dxVec using dynamic field names
                g.vts.dxVec(startN:startN+g.vts.fsdn{ndx,4}-1) = g.(f1).(f2)(1:g.vts.fsdn{ndx,4}, k);
                startN = startN + g.vts.fsdn{ndx,4}; % increment starting dxVec index
            end
        else
            % probably not executed...
            disp(['skpping... ', g.vts.fsdn(ndx,1), g.vts.fsdn(ndx,3) ] );
        end
    end
    
    %dxVec % debug output
end

if flag == 2
    %% Write slnVec vector values to states at index k
    % input =   slnVec - internally generated by the ODE sovler
    %           k - index to write collected states to
    % output = none
    
    %fprintf('DEBUG: Performing g.x.state(:,j) = g.x.dxstate(:,k)...\n')
    
    startN = 1;
    for ndx=1:size(g.vts.fsdn,1)
        % If the number of allocated states is larger than zero
        if g.vts.fsdn{ndx, 4} ~= 0
            % collect field and state name
            f1 = g.vts.fsdn{ndx,1};
            f2 = g.vts.fsdn{ndx, 2};
            
            if strcmp(f1, 'agc')
                for n = 1:g.agc.n_agc
                    % write new states to sace
                    g.agc.agc(n).sace(k) = slnVec(startN);
                    startN = startN+1;
                    for cg = 1:g.agc.agc(n).n_ctrlGen
                        % write new states to controlled gen x
                        g.agc.agc(n).ctrlGen(cg).x(k)  = slnVec(startN);
                        startN = startN +1;
                    end
                end
                
            elseif strcmp(f2, 'pwrmod_p_sigst') || strcmp(f2, 'pwrmod_q_sigst') || ... % pwrmod handling of variable cell states
                   strcmp(f2, 'ivmmod_d_sigst') || strcmp(f2, 'ivmmod_e_sigst') % ivmmod handling of variable cell states
                for n = 1:g.vts.fsdn{ndx,4}
                    % write solution vector to states
                    g.(f1).(f2){n}(k) = slnVec(startN);
                    startN = startN + 1; % increment starting slnVec index
                end
                
            else
                
                for i=1:g.vts.fsdn{ndx, 4}
                    g.(f1).(f2)(i,k) = slnVec(startN+i-1); % write to state col k
                end
                startN = startN+g.vts.fsdn{ndx, 4};
            end
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
            
            
            if strcmp(f1, 'agc')
                for n = 1:g.agc.n_agc
                    % write new states to sace
                    g.vts.stVec(startN) = g.agc.agc(n).sace(k);
                    startN = startN+1;
                    for cg = 1:g.agc.agc(n).n_ctrlGen
                        % write new states to controlled gen x
                        g.vts.stVec(startN) = g.agc.agc(n).ctrlGen(cg).x(k);
                        startN = startN +1;
                    end
                end
                
            elseif strcmp(f2, 'pwrmod_p_sigst') || strcmp(f2, 'pwrmod_q_sigst') || ... % pwrmod handling of variable cell states
                   strcmp(f2, 'ivmmod_d_sigst') || strcmp(f2, 'ivmmod_e_sigst') % ivmmod handling of variable cell states
                for n = 1:g.vts.fsdn{ndx,4}
                    % place states into stVec using dynamic field names
                    g.vts.stVec(startN) = g.(f1).(f2){n}(k);
                    startN = startN + 1; % increment starting stVec index
                end
                
            else
                
                % place state into stVec using dynamic field names
                g.vts.stVec(startN:startN+g.vts.fsdn{ndx,4}-1) = g.(f1).(f2)(1:g.vts.fsdn{ndx,4}, k);
                startN = startN + g.vts.fsdn{ndx,4}; % increment starting stVec index
            end
        else
            % probably not executed...
            disp(['skpping... ', g.vts.fsdn(ndx,1), g.vts.fsdn(ndx,3) ] );
        end
    end
    
end

if flag == 4
    %% DEBUG verify write
    % input ...
    % include at all? - Not updated to accomodate for agc...
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