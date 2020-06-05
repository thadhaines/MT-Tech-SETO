%% Script to detect 'legacy' data input and adjust to new global g approach
% Usage of script to more easily detect legacy global definitons
% global g is already defined
%
%   History:
%   Date        Time    Engineer        Description
%   06/02/20    08:51   Thad Haines     init - lmod_con
%   06/05/20    09:53   Thad Haines     addition of tg_con

% lmod
if exist('lmod_con','var')
    g.lmod.lmod_con = lmod_con;
    clear lmod_con 
end
if ~isfield('g.lmod','lmod_con')
    g.lmod.lmod_con = [];
end

% tg
if exist('tg_con','var')
    % legacy defined
    g.tg.tg_con = tg_con;
    clear tg_con 
end
if ~isfield('g.tg','tg_con')
    % create empty tg_con if not defined
    g.tg.tg_con = [];
end