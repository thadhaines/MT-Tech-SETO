function lm_indx()
% LM_INDX determines the relationship betweel lmod and nc loads.
% LM_INDX determines the relationship between lmod and nc loads by 
% checking for lmod_con and then determining the number of modulated loads
%
%   Syntax:
%   lm_indx()
%
%   History:
%   Date        Time    Engineer        Description
%   08/15/97    17:02   Graham Rogers   Version 1
%   (c) Copyright 1991-1997 Joe H. Chow/ Cherry Tree Scientific Software - All Rights Reserved
%   06/15/20    15:58   Thad Haines     Revised format of globals and internal function documentation

global g % thad
global load_con % required as there is a check for non-conforming loads

% set initial number of load modulations to zero
g.lmod.n_lmod = 0;
% initialize empty index array of loads to modulate
g.lmod.lmod_idx = [];
% handle no modulation case
if ~isfield(g.lmod, 'lmod_con')
    fprintf('*** No load modulation definitions detected.\n')
    g.lmod.lmod_con =[];
end
% checks global lmod_con for load modulation definitions
if ~isempty(g.lmod.lmod_con)
    % collect number of load to modulate
    g.lmod.n_lmod = length(g.lmod.lmod_con(:,1));
    % initialize index array for modulated loads
    g.lmod.lmod_idx = zeros(g.lmod.n_lmod,1);
    
    for j = 1:g.lmod.n_lmod % for each modulated load...
        % find index in load_con
        index = find(g.lmod.lmod_con(j,2) == load_con(:,1));
        
        % if index is found, save index in lmod_idx
        if ~isempty(index)
            g.lmod.lmod_idx(j) = index;
        % else, throw error
        else
            error('*** Load modulation bus must be declared as a non-conforming load.')
        end % end ~isempty(index)
        
    end % end for loop of n_lmod
    
end % end ~isempty(g.lmod.lmod_con)

