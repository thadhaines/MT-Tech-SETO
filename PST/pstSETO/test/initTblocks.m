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

global g

g.vts.t_block = zeros(size(g.sys.sw_con,1)-1, 2);

for n = 1:size(g.sys.sw_con,1)-1
    g.vts.t_block(n,1) = g.sys.sw_con(n,1); % start time
    g.vts.t_block(n,2) = g.sys.sw_con(n+1,1) - g.sys.sw_con(n+1,7); % end time minus timestep
end

g.vts.t_block(end,2) = g.sys.sw_con(end,1); % ensure simulation end time correct
g.vts.t_blockN = []; % init time block index

end