function rlm_indx()
% RLM_INDX determines the relationship betweel rlmod and nc loads.
% RLM_INDX determines the relationship between rlmod and nc loads by
% checking for rlmod_con and then determining the number of modulated loads
%
%   Syntax:
%   rlm_indx()
%
%   History:
%   Date        Time    Engineer        Description
%   08/27/97    17:28   Graham Rogers   Version 1
%   (c) Copyright 1991-1997 Joe H. Chow/ Cherry Tree Scientific Software - All Rights Reserved
%   06/15/20    15:58   Thad Haines     Revised format of globals and internal function documentation
%   07/02/20    13:08   Thad Haines     complete conversion to global g

global g

g.rlmod.n_rlmod = 0;
g.rlmod.rlmod_idx = [];

if ~isempty(g.rlmod.rlmod_con)
    g.rlmod.n_rlmod = length(g.rlmod.rlmod_con(:,1));
    g.rlmod.rlmod_idx = zeros(g.rlmod.n_rlmod,1);
    for j = 1:g.rlmod.n_rlmod
        index = find(g.rlmod.rlmod_con(j,2) == g.ncl.load_con(:,1));
        if ~isempty(index)
            g.rlmod.rlmod_idx(j) = index;
        else
            error('*** The reactive load modulation bus must be declared as a non-conforming load.')
        end
    end
end

