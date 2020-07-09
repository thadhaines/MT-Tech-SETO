function tcsc_indx()
%TCSC_INDX determines the relationship between tcsc and nc loads
% TCSC_INDX determines the relationship between tcsc and nc loads,
% checks for tcsc, dertermines number of tcscs and checks for
% user defined damping controls.
%
% Syntax: tcsc_indx()
%
%   Input: 
%   VOID
%
%   Output: 
%   VOID
%
%   History:
%   Date        Time    Engineer        Description
%   12/29/98    10:29   Graham Rogers   Version 1.0
%   07/09/20    12:03   Thad Haines     Revised format of globals and internal function documentation

global g

if ~isempty(g.tcsc.tcsc_con)
   g.tcsc.n_tcsc = size(g.tcsc.tcsc_con,1);
   %tcsc_idx = zeros(n_tcsc,1);
   for j = 1:g.tcsc.n_tcsc
      index = find(g.tcsc.tcsc_con(j,2)==g.ncl.load_con(:,1));
      if ~isempty(index)
         g.tcsc.tcscf_idx(j) = index;
      else
         error('you must have the tcsc buses declared as a non-conforming load')
      end
      index = find(g.tcsc.tcsc_con(j,3)== g.ncl.load_con(:,1));
      if ~isempty(index)
         g.tcsc.tcsct_idx(j) = index;
      else
         error('you must have the tcsc buses declared as a non-conforming load')
      end
   end
   %check for user defined controls
   if ~isempty(g.tcsc.tcsc_dc)
      g.tcsc.n_tcscud = size(g.tcsc.tcsc_dc,1);
      for j = 1:g.tcsc.n_tcscud
         g.tcsc.dtcscud_idx(j) = g.tcsc.tcsc_dc{j,2};
      end
   end      
end

