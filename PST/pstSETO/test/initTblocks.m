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
%   07/28/20    08:55   Thad Haines     Version 1.0.1 - Allowed time block overlap
%%
global g

g.vts.t_block = zeros(size(g.sys.sw_con,1)-1, 2);
g.vts.fts = {zeros(size(g.sys.sw_con,1)-1, 2)}; % for holding any fixed time blocks

for n = 1:size(g.sys.sw_con,1)-1
    g.vts.t_block(n,1) = g.sys.sw_con(n,1); % start time
    g.vts.t_block(n,2) = g.sys.sw_con(n+1,1);%
    
    % handle optional fixed time blocks
    if ~isempty(g.vts.solver_con)
        if strcmp(g.vts.solver_con{n}, 'huens')
            % remove last time step between blocks for huens to huens
            g.vts.fts{n} = g.vts.t_block(n,1):g.sys.sw_con(n,7):g.vts.t_block(n,2)-g.sys.sw_con(n,7);
            
%             % if another time block exists past the current time blcok
%             if n+1 <= size(g.sys.sw_con,1)-1
%                 % if next time block NOT huens, do not remove last time step
%                 if ~strcmp(g.vts.solver_con{n+1}, 'huens')
%                     g.vts.fts{n} = g.vts.t_block(n,1):g.sys.sw_con(n,7):g.vts.t_block(n,2);
%                 end
%             end
            
        else
            g.vts.fts{n} = 0;
        end
    end
    
end

g.vts.t_block(end,2) = g.sys.sw_con(end,1); % ensure simulation end time correct for variable sims
g.vts.fts{end}(end) = g.sys.sw_con(end,1); % ensure simulation end time correct for fixed step sims
g.vts.t_blockN = []; % init time block index

end