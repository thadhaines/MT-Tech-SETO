function dc_indx(bus,line,dci_dc,dcr_dc)
%DC_INDX forms indexes for the rectifier and inverter in dc load flow.
% DC_INDX  forms indexes for the rectifier and inverter in the dc load flow
%          and indicates the ac buses contected to the converters.
%
% Syntax: dc_indx(bus,line,dci_dc,dcr_dc)
%
%   NOTES:  dci_dc and dcr_dc are the same as g.dc.dci_dc g.dcr_dc
% 
%   Input: 
%   bus - solved loadflow bus data
%   line - line data
%   dcr_dc - user defined damping control at rectifier cell
%   dci_dc - user defined damping control at inverter cell
%
%   Output: 
%   VOID
%
%   History:
%   Date        Time    Engineer        Description
%   01/xx/96    xx:xx   Graham Rogers   Version 1.0
%   01/xx/99    XX:XX   Graham Rogers  	VAddition of user defined damping controls
%   (c) Copyright Joe Chow 1996 - All right reserved
%   07/15/20    10:48   Thad Haines     Revised format of globals and internal function documentation

global g
% pick out ac voltages (should be LT converter transformer buses)
% check that line and cont data is consistent

g.dc.n_conv = 0;
g.dc.n_dcl = 0;
g.dc.dcrud_idx = [];
g.dc.dciud_idx = [];
g.dc.ndcr_ud = 0;
g.dc.ndci_ud = 0;

if ~isempty(g.dc.dcsp_con)
   lconv = length(g.dc.dcsp_con(:,1));
   lline = length(g.dc.dcl_con(:,1));
   lcon = length(g.dc.dcc_con(:,1));
   if lcon~=2*lline || lcon~=lconv % replaced |
      nc = num2str(lconv);
      nl = num2str(lline);
      ncon = num2str(lcon); 
      disp('number of converters = ',nc)
      disp('number of lines = ',nl)
      disp('number of controls = ',ncon)
      error('dc converter and line data inconsistent')
   end
   % find index of rectifier buses
   g.dc.r_idx = find(g.dc.dcsp_con(:,3)==1);
   % find index of inverter buses
   g.dc.i_idx = find(g.dc.dcsp_con(:,3)==2);
   % find index of recitifier current control
   g.dc.ric_idx = find(g.dc.dcc_con(:,9)==1);
   % find index of rectifier power control
   g.dc.rpc_idx = find(g.dc.dcc_con(:,9)==2);
   g.dc.n_dcl = lline;
   g.dc.n_conv = lconv;
   g.dc.inv_ac_bus = g.bus.bus_int(g.dc.dcsp_con(g.dc.i_idx,2));
   g.dc.rec_ac_bus = g.bus.bus_int(g.dc.dcsp_con(g.dc.r_idx,2));
   g.dc.ac_bus = g.bus.bus_int(g.dc.dcsp_con(:,2));
   g.dc.inv_ac_line = zeros(lline,1);
   g.dc.rec_ac_line = zeros(lline,1);
   for j = 1:lline
      acilj = find(g.bus.bus_int(line(:,2)) == g.dc.inv_ac_bus(j));
      if isempty(acilj)
         error(' the inverter bus is not declared as a to bus')
      else
         g.dc.inv_ac_line(j) = acilj;
         acilj = [];
      end
      acrlj = find(g.bus.bus_int(line(:,2)) == g.dc.rec_ac_bus(j));
      if isempty(acrlj)
         error(' the rectifier bus is not declared as a to bus')
      else
         g.dc.rec_ac_line(j) = acrlj;
         acrlj = [];
      end
   end
   g.dc.ac_line = [g.dc.rec_ac_line;g.dc.inv_ac_line];
   % form index of dc lines associated with the inverters
   g.dc.dcli_idx = zeros(g.dc.n_dcl,1);
   for k = 1:g.dc.n_dcl
      g.dc.dcli_idx = g.dc.dcli_idx | (g.dc.dcl_con(k,2)==g.dc.dcsp_con(g.dc.i_idx,1));
   end
   g.dc.dcli_idx = find(g.dc.dcli_idx~=0);
end
g.dc.no_cap_idx = find(g.dc.dcl_con(:,5)==0);
g.dc.cap_idx = find(g.dc.dcl_con(:,5)~=0);
g.dc.l_no_cap = 0;
if ~isempty(g.dc.no_cap_idx); 
   g.dc.l_no_cap = length(g.dc.no_cap_idx);
end
g.dc.l_cap = g.dc.n_dcl-g.dc.l_no_cap;
g.dc.no_ind_idx = find(g.dc.dcl_con(:,4) ==0|g.dc.dcl_con(:,6)==0|g.dc.dcl_con(:,7)==0);

% index of converters in load_con
j = g.bus.bus_int(g.ncl.load_con(:,1));
for k = 1: g.dc.n_conv
   g.dc.ldc_idx(k) = find(j==g.dc.ac_bus(k));
   if isempty(g.dc.ldc_idx(k))
      error('dc converter LT buses must be declared in load_con')
      % an additional load is not allowed at a converter bus
   end
end
% j(ldc_idx(k)) gives  the internal bus number of the kth converter 
% check for user defined controls
[g.dc.ndcr_ud,dummy] = size(g.dc.dcr_dc);
for j = 1:g.dc.ndcr_ud
    if ~isempty(g.dc.dcr_dc{j})
        g.dc.dcrud_idx(j) = g.dc.dcr_dc{j,2};
    end
    g.dc.dcrud_idx(j) = 0;
end     
[g.dc.ndci_ud,dummy] = size(g.dc.dci_dc);
for j = 1:g.dc.ndci_ud
    if ~isempty(g.dc.dci_dc{j})
      g.dc.dciud_idx(j) = g.dc.dci_dc{j,2};
    end
  g.dc.dciud_idx(j) = 0;
end     

return