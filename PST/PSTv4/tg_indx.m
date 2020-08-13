function tg_indx()
%TG_INDX function to create turbine govenor indexes
%
%   Syntax:
%   tg_indx()
%
%   Inputs:
%   NONE
%
%   Output:
%   VOID
%
%   History:
%   Date        Time    Engineer        Description
%   08/15/97    16:58   Graham Rogers   Version 1
%   (c) Copyright 1991-1997 Joe H. Chow/ Cherry Tree Scientific Software - All Rights Reserved
%   06/05/20    09:49   Thad Haines     Revised format of globals and internal function documentation

global g

% initalize with no governos
g.tg.n_tg=0;
g.tg.n_tgh=0;

if ~isempty(g.tg.tg_con)
    
    % Find number of model 1 governors
    g.tg.tg_idx = find(g.tg.tg_con(:,1)==1);
    if ~isempty(g.tg.tg_idx)
        g.tg.n_tg = length(g.tg.tg_idx);
    end
    
    % Find number of model 2 governors (hydro)
    g.tg.tgh_idx = find(g.tg.tg_con(:,1)==2);
    if ~isempty(g.tg.tgh_idx)
        g.tg.n_tgh = length(g.tg.tgh_idx);
    end
end