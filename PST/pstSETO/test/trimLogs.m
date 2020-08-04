function trimLogs(k)
%TRIMLOGS removes excess data allocated for variable time step simulation.
% TRIMLOGS trims logged data to input index k.
%
% Syntax: trimLogs(k)
%
%   NOTES:
%
%   Input:
%   k - data index
%
%   Output:
%   VOID
%
%   History:
%   Date        Time    Engineer        Description
%   07/28/20    11:47   Thad Haines     Version 1
%   08/03/20    15:30   Thad Haines     V 1.0.1 - addition of AGC, area values, and tg_sig

global g
%%
nCell = { ...
    % field, logged values
    'svc', {'B_cv', 'dB_cv', 'svc_sig', 'svc_dsig', 'B_con', 'dB_con', 'xsvc_dc', 'dxsvc_dc' };
    'tcsc', {'B_tcsc', 'dB_tcsc', 'tcsc_sig', 'tcsc_dsig',  'xtcsc_dc', 'dxtcsc_dc', 'td_sig'};
    % '', {}; dc skipped....
    'bus', {'bus_v', 'theta'};
    'lmod', {'lmod_st', 'dlmod_st', 'lmod_sig'};
    'rlmod', {'rlmod_st', 'drlmod_st', 'rlmod_sig'};
    'tg', {'tg1', 'tg2', 'tg3', 'tg4', 'tg5', 'dtg1', 'dtg2', 'dtg3', 'dtg4', 'dtg5', 'tg_sig'};
    'exc', {'V_B', 'exc_sig', 'V_TR', 'V_R', 'V_A', 'V_As', 'Efd', 'R_f', 'dV_TR', 'dV_R', 'dV_As', 'dEfd', 'dR_f'};
    'mac', {'cur_re',   'cur_im',   'psi_re',   'psi_im',  'mac_ang',  'mac_spd', 'dmac_ang', 'dmac_spd',    'pmech',   'pelect',  'edprime',  'eqprime', 'dedprime', 'deqprime',    'psikd',    'psikq',   'dpsikd',   'dpsikq',   'pm_sig',     'curd',     'curq',    'curdg',    'curqg',   'fldcur',       'ed',       'eq',    'eterm',   'qelect',      'vex'};
    'pwr', {'pwrmod_p_st', 'dpwrmod_p_st', 'pwrmod_p_sig',  'pwrmod_q_st', 'dpwrmod_q_st', 'pwrmod_q_sig'};
    'sys', {'t', 'aveF', 'totH'};
    'pss', {'pss1', 'pss2', 'pss3', 'dpss1', 'dpss2', 'dpss3', 'pss_out'};
    'igen', {'vdpig', 'vqpig', 'slig', 'dvdpig', 'dvqpig', 'dslig', 's_igen', 'pig', 'qig', 'tmig'};
    'ind', {'vdp', 'vqp', 'slip', 'dvdp', 'dvqp', 'dslip', 's_mot', 'p_mot', 'q_mot' };
    'vts', {'slns'};
    'area', {'totH', 'aveF', 'totGen', 'totLoad', 'icA', 'icS'}; % add optionally?
    'lmon', {'iFrom', 'iTo', 'sFrom', 'sTo'}; % add optionally?
    'agc', {'Handled Differencly'} ; % add optionally?
    %'', {};
    };

for f=1:size(nCell,1)
    
    % Handle AGC structure
    if strcmp(nCell{f,1}, 'agc')
        aceVars = {'race', 'd_sace', 'sace', 'ace2dist', 'aceSig'};
        cgVars = {'input', 'dx' ,'x'};
        for n = 1:g.agc.n_agc
            % Trim AGC variables
            for sf =1:max(size(aceVars))
                % for each sub-field
                g.agc.agc(n).(aceVars{sf}) = g.agc.agc(n).(aceVars{sf})(:, 1: k);
            end
            
            for cg = 1:g.agc.agc(n).n_ctrlGen
                % Trim CTRL gen variables
                for sf =1:max(size(cgVars))
                    % for each sub-field
                    g.agc.agc(n).ctrlGen(cg).(cgVars{sf}) = g.agc.agc(n).ctrlGen(cg).(cgVars{sf})(:, 1: k);
                end
                
                
            end
        end
    elseif strcmp(nCell{f,1}, 'area')
        % handle logged area values
        for n=1:g.area.n_area
            for sf =1:max(size(nCell{f,2}))
                % for each sub-field
                g.area.area(n).(nCell{f,2}{sf}) = g.area.area(n).(nCell{f,2}{sf})(:, 1: k);
            end
        end
    elseif strcmp(nCell{f,1}, 'lmon')
        % handle logged line monitor values
        for n=1:g.lmon.n_lmon
            for sf =1:max(size(nCell{f,2}))
                % for each sub-field
                g.lmon.line(n).(nCell{f,2}{sf}) = g.lmon.line(n).(nCell{f,2}{sf})(:, 1: k);
            end
        end
    else
        % normal global operation
        % for each row (field) in the nCell...
        for sf =1:max(size(nCell{f,2}))
            % for each sub-field
            g.(nCell{f,1}).(nCell{f,2}{sf}) = g.(nCell{f,1}).(nCell{f,2}{sf})(:, 1: k);
        end
    end
end

end