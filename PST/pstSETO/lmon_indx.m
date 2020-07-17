function lmon_indx
%LMON_INDX checks if line monitoring indicies exist
% LMON_INDX checks if line monitoring indicies exist
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
    if all(g.lmon.lmon_con <= size(g.line.lineOG,1))
        g.lmon.n_lmon = length(g.lmon.lmon_con); % count lines to monitor
        g.lmon.busFromTo = zeros(g.lmon.n_lmon,2); % array for collecting bus info
        
        for Lndx = 1:g.lmon.n_lmon
           % create placeholder structure for line flow info
           g.lmon.line(Lndx).busArrayNdx = g.lmon.n_lmon(Lndx);
           g.lmon.line(Lndx).FromBus = g.line.lineOG(g.lmon.lmon_con(Lndx), 1);
           g.lmon.line(Lndx).ToBus = g.line.lineOG(g.lmon.lmon_con(Lndx), 2);
           
           % vector of busses for later calculations
           g.lmon.busFromTo(Lndx,1) = g.line.lineOG(g.lmon.lmon_con(Lndx), 1);
           g.lmon.busFromTo(Lndx,2) = g.line.lineOG(g.lmon.lmon_con(Lndx), 2);
           
           % placeholders for logged data
           g.lmon.line(Lndx).iFrom = [];
           g.lmon.line(Lndx).iTo = [];
           g.lmon.line(Lndx).sFrom = [];
           g.lmon.line(Lndx).sTo = [];
        end
    else
        error('Line monitoring index outside of line array')
    end
else
    g.lmon.n_lmon = 0;
end