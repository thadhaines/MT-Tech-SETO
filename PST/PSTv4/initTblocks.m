function initTblocks
%INITTblocks creates time blocks for vts
% INITSTEP creates time blocks for vts from global g.sys.sw_con
%
% Syntax: initTblocks
%
%   NOTES:
%
%   Input:
%   VOID
%
%   Output:
%   VOID
%
%   History:
%   Date        Time    Engineer        Description
%   07/27/20    11:06   Thad Haines     Version 1
%   07/28/20    08:55   Thad Haines     Version 1.0.1 - fixed time block overlap
%   08/14/20    09:47   Thad Haines     Version 1.0.2 - fixed time step now handled same as PST3
%%
global g

g.vts.t_block = zeros(size(g.sys.sw_con,1)-1, 2);
g.vts.fts = {zeros(size(g.sys.sw_con,1)-1, 2)}; % for holding any fixed time blocks
g.vts.fts_dc = {zeros(size(g.sys.sw_con,1)-1, 2)}; % for holding any fixed DC time blocks

for n = 1:size(g.sys.sw_con,1)-1
    g.vts.t_block(n,1) = g.sys.sw_con(n,1); % start time
    g.vts.t_block(n,2) = g.sys.sw_con(n+1,1);% end time
    
    % handle optional fixed time blocks
    if ~isempty(g.vts.solver_con)
        if strcmp(g.vts.solver_con{n}, 'huens')
            % mimic original PST time vector creation
            
            if g.sys.sw_con(n,7) == 0
                g.sys.sw_con(n,7) = 0.01;   % enusre time step not 0
            end
            
            nSteps = fix((g.vts.t_block(n,2)-g.vts.t_block(n,1))/g.sys.sw_con(n,7));
            
            if nSteps == 0
                nSteps = 1; % ensure at least 1 step is taken
            end
            
            g.k.h(n) = (g.vts.t_block(n,2)-g.vts.t_block(n,1))/nSteps; % adjust timestep
            g.k.h_dc(n) = g.k.h(n)/10; % h_dc 10 times faster
            
            if n+1 <= size(g.vts.solver_con,1)              % if another time block exists,
                if strcmp(g.vts.solver_con{n+1}, 'huens')   % and is huens method
                    % remove last time step between blocks for huens to huens
                    g.vts.fts{n} = g.vts.t_block(n,1):g.k.h(n):g.vts.t_block(n,2)-g.k.h(n);
                else
                    % Create time block with last time step
                    g.vts.fts{n} = g.vts.t_block(n,1):g.k.h(n):g.vts.t_block(n,2);
                end
            else
                % Create time block with last time step
                g.vts.fts{n} = g.vts.t_block(n,1):g.k.h(n):g.vts.t_block(n,2);
            end
            
            % DC time vector...
            if ~isempty(g.dc.dcl_con)
                g.vts.fts_dc{n} = g.vts.t_block(n,1):g.k.h_dc(n):g.vts.t_block(n,2)-g.k.h_dc(n);
            else
                g.vts.fts_dc{n} = 0;
            end
            
        else
            g.vts.fts{n} = 0;
            g.vts.fts_dc{n} = 0;
        end
    end
    
end

%g.vts.t_block(end,2) = g.sys.sw_con(end,1); % ensure simulation end time correct for variable sims
%g.vts.fts{end}(end) = g.sys.sw_con(end,1); % ensure simulation end time correct for fixed step sims
g.vts.t_blockN = []; % init time block index

end