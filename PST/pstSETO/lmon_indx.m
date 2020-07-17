function lmon_indx
%LMON_INDX creates required indicies for line monitoring
% LMON_INDX creates required indicies for line monitoring
%
% Syntax: lmon_indx
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
%   07/17/20    07:59   Thad Haines     Version 1

global g

if ~isempty(g.lmon.lmon_con)
    % ensure lines are available in line_con to monitor
   g.lmon.n_lmon = size(g.lmon.lmon_con);
   
else
    g.lmon.n_lmon = 0;
end