function lm_indx()
% syntax: f = lm_indx
% 5:02 PM 15/08/97
% determines the relationship between lmod and nc loads
% checks for lmod
% determines number of modulated loads

% original PST code
% f is a dummy variable
% f = 0;
% global lmod_con load_con  n_lmod  lmod_idx
% n_lmod = 0;
% lmod_idx = [];
% if ~isempty(lmod_con)
%     n_lmod = length(lmod_con(:,1));
%     lmod_idx = zeros(n_lmod,1);
%     for j = 1:n_lmod
%        index = find(lmod_con(j,2)==load_con(:,1));
%        if ~isempty(index)
%           lmod_idx(j) = index;
%        else
%           error('you must have the load modulation bus declared as a non-conforming load')
%        end
%     end
% end

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

