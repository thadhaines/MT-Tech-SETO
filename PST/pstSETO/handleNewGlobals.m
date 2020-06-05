%% Script to detect 'legacy' input and adjust to new global g approach
% global g is already defined
%
%   History:
%   Date        Time    Engineer        Description
%   06/02/20    08:51   Thad Haines     init - lmod

% lmod
if exist('lmod_con','var')
    g.lmod.lmod_con = lmod_con;
    clear lmod_con 
end